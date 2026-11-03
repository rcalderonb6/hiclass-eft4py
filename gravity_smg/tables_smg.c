/** @file tables_smg.c Tabulated inputs for hi_class parametrisations
 *
 * Reading, validation, interpolation and extrapolation of
 *  - gravity_model = alphas_table (alpha_K, alpha_B, alpha_M, alpha_T, Delta M_*^2 vs ln a)
 *  - expansion_model = table      (rho_smg/rho_smg(a=1) [, w_smg] vs ln a)
 *
 * See include/tables_smg.h for the conventions.
 */

#include "tables_smg.h"
#include "input.h"
#include "background.h"

#include <errno.h>
#include <math.h>

/* Tolerance on ln a when a table is evaluated just past its last node
   (hi_class evaluates exactly at a = 1, up to round-off). */
#define _TABLE_LNA_EPS_SMG_ 1e-10
/* Early-time continuation: tabulated down to ln a = ln(1e-30), step 0.05 */
#define _TABLE_LNA_LOW_SMG_ (-69.07755278982137)
#define _TABLE_EXT_STEP_SMG_ 0.05

/* ------------------------------------------------------------------------ */
/* Reading                                                                  */
/* ------------------------------------------------------------------------ */

/**
 * Read a table from a file or from comma-separated lists.
 *
 * @param pfc        Input: file content
 * @param file_key   Input: name of the file key
 * @param list_keys  Input: names of the list keys, one per column (column 0 is ln a)
 * @param n_cols     Input: total number of columns, including ln a
 * @param n_required Input: number of leading columns that must be present
 * @param n          Output: number of nodes
 * @param n_found    Output: number of columns found (n_required <= n_found <= n_cols)
 * @param cols       Output: n_cols arrays (NULL for missing columns), allocated here
 * @param errmsg     Output: error message
 */
static int tables_smg_read_columns(struct file_content * pfc,
                                   char * file_key,
                                   char ** list_keys,
                                   int n_cols,
                                   int n_required,
                                   int * n,
                                   int * n_found,
                                   double ** cols,
                                   ErrorMsg errmsg) {

  int j, flag_file, flag_list, any_list, size;
  FileArg * file_name;
  double * list;

  for (j = 0; j < n_cols; j++) cols[j] = NULL;
  *n = 0;
  *n_found = 0;

  /* FileArg is large (_ARGUMENT_LENGTH_MAX_): keep it off the stack */
  class_alloc(file_name, sizeof(FileArg), errmsg);
  class_call(parser_read_string(pfc, file_key, file_name, &flag_file, errmsg),
             errmsg, errmsg);

  /* lists */
  any_list = _FALSE_;
  for (j = 0; j < n_cols; j++) {
    class_call(parser_read_list_of_doubles(pfc, list_keys[j], &size, &list, &flag_list, errmsg),
               errmsg, errmsg);
    if (flag_list == _TRUE_) {
      any_list = _TRUE_;
      cols[j] = list;
      if (j == 0) {
        *n = size;
      }
      else {
        class_test(cols[0] == NULL, errmsg,
                   "'%s' given without '%s'", list_keys[j], list_keys[0]);
        class_test(size != *n, errmsg,
                   "'%s' has %d entries, but '%s' has %d: all table columns must have the same length",
                   list_keys[j], size, list_keys[0], *n);
      }
    }
  }

  class_test((flag_file == _TRUE_) && (any_list == _TRUE_), errmsg,
             "give the table either as a file ('%s') or as lists ('%s', ...), not both",
             file_key, list_keys[0]);
  class_test((flag_file == _FALSE_) && (any_list == _FALSE_), errmsg,
             "no table given: set '%s', or the lists '%s', ...", file_key, list_keys[0]);

  if (any_list == _TRUE_) {
    for (j = 0; j < n_required; j++)
      class_test(cols[j] == NULL, errmsg, "table column '%s' is required", list_keys[j]);
    *n_found = n_required;
    for (j = n_required; j < n_cols; j++) {
      if (cols[j] != NULL) {
        class_test(*n_found != j, errmsg,
                   "optional table column '%s' given without the preceding ones", list_keys[j]);
        (*n_found)++;
      }
    }
    free(file_name);
    return _SUCCESS_;
  }

  /* file: one row per line, blank-separated, '#' starts a comment */
  {
    FILE * f;
    char line[_LINE_LENGTH_MAX_];
    int capacity = 256, row_cols, line_number = 0;
    char * cursor, * end;
    double value;

    f = fopen(*file_name, "r");
    class_test(f == NULL, errmsg, "could not open table file '%s' (%s)", *file_name, strerror(errno));

    for (j = 0; j < n_cols; j++)
      class_alloc(cols[j], capacity*sizeof(double), errmsg);

    while (fgets(line, _LINE_LENGTH_MAX_, f) != NULL) {
      line_number++;
      class_test(strchr(line, '\n') == NULL && !feof(f), errmsg,
                 "line %d of '%s' is longer than %d characters", line_number, *file_name, _LINE_LENGTH_MAX_-1);
      if ((cursor = strchr(line, '#')) != NULL) *cursor = '\0';

      row_cols = 0;
      cursor = line;
      while (1) {
        while (*cursor == ' ' || *cursor == '\t' || *cursor == ',' || *cursor == '\r' || *cursor == '\n') cursor++;
        if (*cursor == '\0') break;
        errno = 0;
        value = strtod(cursor, &end);
        class_test(end == cursor, errmsg,
                   "line %d of '%s': cannot read a number at '%.20s'", line_number, *file_name, cursor);
        class_test(row_cols >= n_cols, errmsg,
                   "line %d of '%s' has more than %d columns", line_number, *file_name, n_cols);
        if (*n == capacity) {
          capacity *= 2;
          for (j = 0; j < n_cols; j++)
            class_realloc(cols[j], capacity*sizeof(double), errmsg);
        }
        cols[row_cols][*n] = value;
        row_cols++;
        cursor = end;
      }
      if (row_cols == 0) continue; /* blank or comment line */

      if (*n_found == 0) {
        class_test(row_cols < n_required, errmsg,
                   "line %d of '%s' has %d columns, at least %d are required",
                   line_number, *file_name, row_cols, n_required);
        *n_found = row_cols;
      }
      class_test(row_cols != *n_found, errmsg,
                 "line %d of '%s' has %d columns, previous lines have %d",
                 line_number, *file_name, row_cols, *n_found);
      (*n)++;
    }
    fclose(f);

    for (j = *n_found; j < n_cols; j++) {
      free(cols[j]);
      cols[j] = NULL;
    }
  }

  free(file_name);
  return _SUCCESS_;
}

/** Read table_extrapolation_smg (default: scaling) and table_extrapolation_window_smg. */
static int tables_smg_read_extrapolation(struct file_content * pfc,
                                         struct table_extrapolation_settings_smg * pset,
                                         ErrorMsg errmsg) {
  FileArg * string;
  int flag;
  double window;

  pset->mode = table_scaling;
  pset->window = log(10.);
  pset->tol = 0.2;
  pset->verbose = 0;
  class_call(parser_read_int(pfc, "background_verbose", &(pset->verbose), &flag, errmsg), errmsg, errmsg);
  if (flag == _FALSE_) pset->verbose = 0;

  class_alloc(string, sizeof(FileArg), errmsg);
  class_call(parser_read_string(pfc, "table_extrapolation_smg", string, &flag, errmsg),
             errmsg, errmsg);
  if (flag == _TRUE_) {
    if (strcmp(*string, "scaling") == 0)
      pset->mode = table_scaling;
    else if (strcmp(*string, "power_law") == 0)
      pset->mode = table_power_law;
    else if (strcmp(*string, "constant") == 0)
      pset->mode = table_constant;
    else
      class_stop(errmsg, "table_extrapolation_smg = '%s': expected 'scaling', 'power_law' or 'constant'", *string);
  }
  free(string);

  class_call(parser_read_double(pfc, "table_extrapolation_tol_smg", &window, &flag, errmsg),
             errmsg, errmsg);
  if (flag == _TRUE_) pset->tol = window;

  class_call(parser_read_double(pfc, "table_extrapolation_window_smg", &window, &flag, errmsg),
             errmsg, errmsg);
  if (flag == _TRUE_) {
    class_test(window <= 0., errmsg, "table_extrapolation_window_smg must be positive, got %g", window);
    pset->window = window;
  }
  return _SUCCESS_;
}

/** Common checks: enough nodes, ln a strictly increasing and reaching a = 1, finite values. */
static int tables_smg_validate(char * what,
                               int n,
                               int n_found,
                               double ** cols,
                               ErrorMsg errmsg) {
  int i, j;

  class_test(n < 4, errmsg, "%s: at least 4 nodes are needed, got %d", what, n);
  for (i = 0; i < n; i++)
    for (j = 0; j < n_found; j++)
      class_test(!isfinite(cols[j][i]), errmsg,
                 "%s: non-finite value in column %d at node %d", what, j, i);
  for (i = 1; i < n; i++)
    class_test(cols[0][i] <= cols[0][i-1], errmsg,
               "%s: ln a must be strictly increasing (node %d: %.17g after %.17g)",
               what, i, cols[0][i], cols[0][i-1]);
  class_test(cols[0][n-1] < -_TABLE_LNA_EPS_SMG_, errmsg,
             "%s: the last node is at ln a = %g < 0; tables must reach a = 1 (no extrapolation to late times)",
             what, cols[0][n-1]);
  return _SUCCESS_;
}

/* ------------------------------------------------------------------------ */
/* Generic table: spline, evaluation, simple tail below the first node      */
/* ------------------------------------------------------------------------ */

/** Two-node tail slope of one column (0 = constant). */
static double tables_smg_tail_slope(double x0, double y0, double x1, double y1,
                                    short is_log, short constant) {
  if (constant == _TRUE_) return 0.;
  if (is_log == _TRUE_) return (y1 - y0)/(x1 - x0);
  if (y0 == 0. || y1 == 0. || (y0 > 0.) != (y1 > 0.)) return 0.;
  return log(y1/y0)/(x1 - x0);
}

/**
 * Fill a struct table_smg from columns (copied), compute the spline and the
 * tail slopes used below the first node.
 */
static int tables_smg_init_table(struct table_smg * ptab,
                                 int n,
                                 int ncol,
                                 double * lna,
                                 double ** cols,
                                 short * is_log,
                                 short constant_tail,
                                 ErrorMsg errmsg) {
  int i, j;

  ptab->n = n;
  ptab->ncol = ncol;
  class_alloc(ptab->lna, n*sizeof(double), errmsg);
  class_alloc(ptab->y, n*ncol*sizeof(double), errmsg);
  class_alloc(ptab->ddy, n*ncol*sizeof(double), errmsg);
  class_alloc(ptab->is_log, ncol*sizeof(short), errmsg);
  class_alloc(ptab->extrap_slope, ncol*sizeof(double), errmsg);
  class_alloc(ptab->sign, ncol*sizeof(short), errmsg);

  for (i = 0; i < n; i++) {
    ptab->lna[i] = lna[i];
    for (j = 0; j < ncol; j++)
      ptab->y[i*ncol+j] = cols[j][i];
  }
  for (j = 0; j < ncol; j++) {
    short nonneg = _TRUE_, nonpos = _TRUE_;
    for (i = 0; i < n; i++) {
      if (cols[j][i] < 0.) nonneg = _FALSE_;
      if (cols[j][i] > 0.) nonpos = _FALSE_;
    }
    ptab->sign[j] = (is_log[j] == _TRUE_) ? 0 : (nonneg == _TRUE_ ? 1 : (nonpos == _TRUE_ ? -1 : 0));
    ptab->is_log[j] = is_log[j];
    ptab->extrap_slope[j] = tables_smg_tail_slope(lna[0], cols[j][0], lna[1], cols[j][1],
                                                  is_log[j], constant_tail);
  }

  class_call(array_spline_table_lines(ptab->lna, n, ptab->y, ncol, ptab->ddy,
                                      _SPLINE_EST_DERIV_, errmsg),
             errmsg, errmsg);
  return _SUCCESS_;
}

static int tables_smg_free_table(struct table_smg * ptab) {
  free(ptab->lna);
  free(ptab->y);
  free(ptab->ddy);
  free(ptab->is_log);
  free(ptab->extrap_slope);
  free(ptab->sign);
  return _SUCCESS_;
}

/**
 * Evaluate all columns of a table at ln a, and optionally their derivatives
 * with respect to ln a (dresult may be NULL).
 */
static int tables_smg_eval_table(struct table_smg * ptab,
                                 double lna,
                                 double * result,
                                 double * dresult,
                                 ErrorMsg errmsg) {
  int j, last_index = 0;
  double dx;

  if (lna < ptab->lna[0]) {
    dx = lna - ptab->lna[0];
    for (j = 0; j < ptab->ncol; j++) {
      double y0 = ptab->y[j], s = ptab->extrap_slope[j];
      if (ptab->is_log[j] == _TRUE_) {
        result[j] = y0 + s*dx;
        if (dresult != NULL) dresult[j] = s;
      }
      else {
        result[j] = y0*exp(s*dx);
        if (dresult != NULL) dresult[j] = s*result[j];
      }
    }
    return _SUCCESS_;
  }

  class_test(lna > ptab->lna[ptab->n-1] + _TABLE_LNA_EPS_SMG_, errmsg,
             "table evaluated at ln a = %g, after its last node ln a = %g", lna, ptab->lna[ptab->n-1]);
  lna = MIN(lna, ptab->lna[ptab->n-1]);

  class_call(array_interpolate_spline(ptab->lna, ptab->n, ptab->y, ptab->ddy, ptab->ncol,
                                      lna, &last_index, result, ptab->ncol, errmsg),
             errmsg, errmsg);
  for (j = 0; j < ptab->ncol; j++) {
    if (ptab->sign[j] > 0) result[j] = MAX(result[j], 0.);
    if (ptab->sign[j] < 0) result[j] = MIN(result[j], 0.);
  }
  if (dresult != NULL) {
    class_call(array_derivate_spline(ptab->lna, ptab->n, ptab->y, ptab->ddy, ptab->ncol,
                                     lna, &last_index, dresult, ptab->ncol, errmsg),
               errmsg, errmsg);
  }
  return _SUCCESS_;
}

/* ------------------------------------------------------------------------ */
/* Early-time continuation: f = f_0 (a/a_0)^s (rho_bg/rho_bg,0)^beta         */
/* ------------------------------------------------------------------------ */

/**
 * ln rho_bg and d ln rho_bg / d ln a, with rho_bg hi_class's radiation +
 * matter density (photons, ur, baryons, cdm, ncdm), in the units of rho_smg.
 */
static int tables_smg_ln_rho_bg(struct background * pba,
                                double lna,
                                double * L,
                                double * dL,
                                ErrorMsg errmsg) {
  double a = exp(lna), H02 = pba->H0*pba->H0;
  double rho_r = H02*(pba->Omega0_g + pba->Omega0_ur)*pow(a,-4);
  double rho_m = H02*(pba->Omega0_b + pba->Omega0_cdm)*pow(a,-3);
  double rho = rho_r + rho_m, drho = -4.*rho_r - 3.*rho_m;
  double rho_ncdm, p_ncdm;
  int n_ncdm;

  for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++) {
    class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                       pba->w_ncdm_bg[n_ncdm],
                                       pba->q_size_ncdm_bg[n_ncdm],
                                       pba->M_ncdm[n_ncdm],
                                       pba->factor_ncdm[n_ncdm],
                                       1./a-1.,
                                       NULL, &rho_ncdm, &p_ncdm, NULL, NULL),
               errmsg, errmsg);
    rho += rho_ncdm;
    drho += -3.*(rho_ncdm + p_ncdm);
  }
  *L = log(rho);
  *dL = drho/rho;
  return _SUCCESS_;
}

/** Number of nodes within the first `window` e-folds (at least 4). */
static int tables_smg_count_window(int n, double * x, double window) {
  int k = 0;
  while ((k < n) && (x[k] <= x[0] + window)) k++;
  return MIN(n, MAX(k, 4));
}

/** Continuation parameters of one column: g(x) = g_0 + s (x - x_0) + beta (L(x) - L_0). */
struct table_scaling_fit_smg {
  short constant; /**< hold the first-node value */
  double s;
  double beta;
  double rms;     /**< rms residual of the fit over the window, in ln f */
};

/**
 * Fit the continuation of one column over the first nodes (pinned to node 0).
 * g = y (is_log) or ln|y| (otherwise). L: ln rho_bg at the nodes.
 */
static struct table_scaling_fit_smg tables_smg_fit_scaling(struct table_extrapolation_settings_smg * pset,
                                                           int n_fit,
                                                           double * x,
                                                           double * L,
                                                           double * y,
                                                           int stride,
                                                           short is_log) {
  struct table_scaling_fit_smg fit = {_FALSE_, 0., 0., 0.};
  double Sxx = 0., SLL = 0., SxL = 0., Sxg = 0., SLg = 0., det;
  double y0 = y[0], g0, g, dx, dL, dg;
  int i;

  if (pset->mode == table_constant) {
    fit.constant = _TRUE_;
    return fit;
  }
  if (is_log == _FALSE_) {
    for (i = 0; i < n_fit; i++) {
      double yi = y[i*stride];
      if (yi == 0. || (yi > 0.) != (y0 > 0.)) {
        fit.constant = _TRUE_;
        return fit;
      }
    }
  }

  g0 = (is_log == _TRUE_) ? y0 : log(fabs(y0));
  for (i = 1; i < n_fit; i++) {
    g = (is_log == _TRUE_) ? y[i*stride] : log(fabs(y[i*stride]));
    dx = x[i] - x[0];
    dL = L[i] - L[0];
    Sxx += dx*dx; SLL += dL*dL; SxL += dx*dL;
    Sxg += dx*(g - g0); SLg += dL*(g - g0);
  }

  det = Sxx*SLL - SxL*SxL;
  if ((pset->mode == table_power_law) || (det <= 1e-10*Sxx*SLL)) {
    /* power law; also when the window is too deep in a single era for
       rho_bg to differ from a power law (s and beta degenerate there) */
    fit.s = Sxg/Sxx;
    fit.beta = 0.;
  }
  else {
    fit.s = (SLL*Sxg - SxL*SLg)/det;
    fit.beta = (Sxx*SLg - SxL*Sxg)/det;
  }
  for (i = 1; i < n_fit; i++) {
    g = (is_log == _TRUE_) ? y[i*stride] : log(fabs(y[i*stride]));
    dg = g - g0 - fit.s*(x[i] - x[0]) - fit.beta*(L[i] - L[0]);
    fit.rms += dg*dg;
  }
  fit.rms = sqrt(fit.rms/(n_fit - 1));
  return fit;
}

/** Value of a continued column at x, given L(x) - L_0. */
static double tables_smg_continued(struct table_scaling_fit_smg * pfit,
                                   double y0, short is_log, double dx, double dL) {
  double dg;
  if (pfit->constant == _TRUE_) return y0;
  dg = pfit->s*dx + pfit->beta*dL;
  return (is_log == _TRUE_) ? y0 + dg : y0*exp(dg);
}

/**
 * Reliability check of one column's continuation: refit over half the window
 * and compare the two continuations at lna_check (see tables_smg.h).
 */
static int tables_smg_check_fit(struct table_extrapolation_settings_smg * pset,
                                char * what,
                                struct table_scaling_fit_smg * pfit,
                                int n, double * x, double * L_fit,
                                double * y, int stride, short is_log,
                                double lna_check, double L_check,
                                ErrorMsg errmsg) {
  struct table_scaling_fit_smg half;
  int n_half;
  double dx = lna_check - x[0], dL = L_check - L_fit[0], delta, sa, ba, sb, bb;

  if ((pset->tol <= 0.) || (pset->mode == table_constant) || (dx >= 0.)) return _SUCCESS_;
  n_half = tables_smg_count_window(n, x, 0.5*pset->window);
  half = tables_smg_fit_scaling(pset, n_half, x, L_fit, y, stride, is_log);
  sa = (pfit->constant == _TRUE_) ? 0. : pfit->s;  ba = (pfit->constant == _TRUE_) ? 0. : pfit->beta;
  sb = (half.constant == _TRUE_) ? 0. : half.s;    bb = (half.constant == _TRUE_) ? 0. : half.beta;
  delta = fabs((sa - sb)*dx + (ba - bb)*dL);
  if (pset->verbose > 1)
    printf(" -> table continuation of %-8s from a = %.3e: s = %+.4f, beta = %+.4f; "
           "half-window difference at a = %.3e: %.2e\n",
           what, exp(x[0]), sa, ba, exp(lna_check), delta);
  class_test(delta > pset->tol, errmsg,
             "the early-time continuation of %s is not determined by the table: fitting over the first "
             "%g or %g e-folds after its first node (a = %.3e) gives values differing by %.2g in ln at "
             "a = %.3e (tolerance table_extrapolation_tol_smg = %g). Start the table earlier, "
             "or loosen the tolerance if this column does not matter at early times.",
             what, pset->window, 0.5*pset->window, exp(x[0]), delta, exp(lna_check), pset->tol);
  return _SUCCESS_;
}

/**
 * Prepare the continuation grid below the first node: number of nodes in
 * the fit window, the continuation nodes x_ext[0..m-1] (< x_0), and
 * L = ln rho_bg at both. Returns m = 0 if no continuation is needed.
 */
static int tables_smg_continuation_grid(struct background * pba,
                                        struct table_extrapolation_settings_smg * pset,
                                        int n,
                                        double * x,
                                        int * n_fit,
                                        double ** L_fit,
                                        int * m,
                                        double ** x_ext,
                                        double ** L_ext,
                                        double ** dL_ext,
                                        ErrorMsg errmsg) {
  int i, k;
  double dummy;

  *n_fit = tables_smg_count_window(n, x, pset->window);

  class_alloc(*L_fit, (*n_fit)*sizeof(double), errmsg);
  for (i = 0; i < *n_fit; i++)
    class_call(tables_smg_ln_rho_bg(pba, x[i], &((*L_fit)[i]), &dummy, errmsg), errmsg, errmsg);

  *m = (int)ceil((x[0] - _TABLE_LNA_LOW_SMG_)/_TABLE_EXT_STEP_SMG_);
  if (*m < 2) *m = 0;
  class_alloc(*x_ext, MAX(*m,1)*sizeof(double), errmsg);
  class_alloc(*L_ext, MAX(*m,1)*sizeof(double), errmsg);
  class_alloc(*dL_ext, MAX(*m,1)*sizeof(double), errmsg);
  for (k = 0; k < *m; k++) {
    (*x_ext)[k] = x[0] - (*m - k)*_TABLE_EXT_STEP_SMG_;
    class_call(tables_smg_ln_rho_bg(pba, (*x_ext)[k], &((*L_ext)[k]), &((*dL_ext)[k]), errmsg), errmsg, errmsg);
  }
  return _SUCCESS_;
}

/* ------------------------------------------------------------------------ */
/* alphas_table                                                             */
/* ------------------------------------------------------------------------ */

int tables_smg_read_alphas(struct file_content * pfc,
                           struct alphas_table_smg * ptab,
                           ErrorMsg errmsg) {

  /* column order in files and lists: ln a, alpha_K, alpha_B, alpha_M, alpha_T, Delta M_*^2 */
  char * keys[6] = {"table_lna_smg", "table_alpha_K_smg", "table_alpha_B_smg",
                    "table_alpha_M_smg", "table_alpha_T_smg", "table_Delta_M2_smg"};
  double * cols[6];
  double * spl_cols[3];
  short is_log[3] = {_FALSE_, _FALSE_, _FALSE_};
  short constant;
  int n, n_found, i, j, flag;
  double tol = 1e-2, scale, fd, mean;

  class_call(tables_smg_read_columns(pfc, "alphas_table_file_smg", keys, 6, 6, &n, &n_found, cols, errmsg),
             errmsg, errmsg);
  class_call(tables_smg_validate("alphas_table", n, n_found, cols, errmsg), errmsg, errmsg);
  class_call(tables_smg_read_extrapolation(pfc, &(ptab->extrapolation), errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "table_alpha_M_tol_smg", &tol, &flag, errmsg), errmsg, errmsg);
  if (flag == _FALSE_) tol = 1e-2;

  /* M_*^2 > 0, and alpha_M consistent with ln M_*^2 between nodes */
  class_alloc(ptab->lnM2, n*sizeof(double), errmsg);
  class_alloc(ptab->alpha_M, n*sizeof(double), errmsg);
  scale = 0.;
  for (i = 0; i < n; i++) {
    class_test(1. + cols[5][i] <= 0., errmsg,
               "alphas_table: M_*^2 = 1 + Delta_M2 = %g <= 0 at node %d", 1. + cols[5][i], i);
    ptab->lnM2[i] = log1p(cols[5][i]);
    ptab->alpha_M[i] = cols[3][i];
    scale = MAX(scale, fabs(cols[3][i]));
  }
  for (i = 1; i < n; i++)
    scale = MAX(scale, fabs((ptab->lnM2[i] - ptab->lnM2[i-1])/(cols[0][i] - cols[0][i-1])));
  for (i = 1; i < n; i++) {
    fd = (ptab->lnM2[i] - ptab->lnM2[i-1])/(cols[0][i] - cols[0][i-1]);
    mean = 0.5*(ptab->alpha_M[i] + ptab->alpha_M[i-1]);
    class_test(fabs(fd - mean) > tol*scale, errmsg,
               "alphas_table: alpha_M is inconsistent with M_*^2 between nodes %d and %d "
               "(d ln M_*^2/d ln a = %g, mean alpha_M = %g, tolerance %g x %g). "
               "Check the columns, or loosen table_alpha_M_tol_smg.",
               i-1, i, fd, mean, tol, scale);
  }

  class_alloc(ptab->dlnM2, (n-1)*sizeof(double), errmsg);
  for (i = 0; i < n-1; i++)
    ptab->dlnM2[i] = ptab->lnM2[i+1] - ptab->lnM2[i];

  constant = (ptab->extrapolation.mode == table_constant) ? _TRUE_ : _FALSE_;
  ptab->alpha_M_slope = tables_smg_tail_slope(cols[0][0], cols[3][0], cols[0][1], cols[3][1],
                                              _FALSE_, constant);

  spl_cols[0] = cols[1];
  spl_cols[1] = cols[2];
  spl_cols[2] = cols[4];
  class_call(tables_smg_init_table(&(ptab->spl), n, 3, cols[0], spl_cols, is_log, constant, errmsg),
             errmsg, errmsg);

  /* the scaling/power-law continuation needs hi_class's background: done at first use */
  ptab->prepared = constant;

  for (j = 0; j < 6; j++) free(cols[j]);
  return _SUCCESS_;
}

/** Add the early-time continuation of the alphas (see tables_smg.h). */
static int tables_smg_prepare_alphas(struct background * pba,
                                     struct alphas_table_smg * ptab,
                                     double lna_check,
                                     ErrorMsg errmsg) {
  struct table_smg * pspl = &(ptab->spl);
  int n = pspl->n, n_fit, m, k, i, j, q, n_new;
  double * L_fit, * x_ext, * L_ext, * dL_ext;
  double * new_x, * new_cols[3], * new_lnM2, * new_aM, * new_dlnM2;
  struct table_scaling_fit_smg fit[3], fit_aM;
  short is_log[3] = {_FALSE_, _FALSE_, _FALSE_};
  /* 3-point Gauss-Legendre on [-1, 1] */
  double gx[3] = {-sqrt(0.6), 0., sqrt(0.6)}, gw[3] = {5./9., 8./9., 5./9.};

  class_call(tables_smg_continuation_grid(pba, &(ptab->extrapolation), n, pspl->lna,
                                          &n_fit, &L_fit, &m, &x_ext, &L_ext, &dL_ext, errmsg),
             errmsg, errmsg);

  if (m > 0) {
    for (j = 0; j < 3; j++)
      fit[j] = tables_smg_fit_scaling(&(ptab->extrapolation), n_fit, pspl->lna, L_fit, pspl->y + j, 3, _FALSE_);
    fit_aM = tables_smg_fit_scaling(&(ptab->extrapolation), n_fit, pspl->lna, L_fit, ptab->alpha_M, 1, _FALSE_);

    {
      char * names[3] = {"alpha_K", "alpha_B", "alpha_T"};
      double L_check, dummy;
      class_call(tables_smg_ln_rho_bg(pba, lna_check, &L_check, &dummy, errmsg), errmsg, errmsg);
      for (j = 0; j < 3; j++)
        class_call(tables_smg_check_fit(&(ptab->extrapolation), names[j], &fit[j], n, pspl->lna, L_fit,
                                        pspl->y + j, 3, _FALSE_, lna_check, L_check, errmsg),
                   errmsg, errmsg);
      class_call(tables_smg_check_fit(&(ptab->extrapolation), "alpha_M", &fit_aM, n, pspl->lna, L_fit,
                                      ptab->alpha_M, 1, _FALSE_, lna_check, L_check, errmsg),
                 errmsg, errmsg);
    }

    n_new = m + n;
    class_alloc(new_x, n_new*sizeof(double), errmsg);
    for (j = 0; j < 3; j++) class_alloc(new_cols[j], n_new*sizeof(double), errmsg);
    class_alloc(new_lnM2, n_new*sizeof(double), errmsg);
    class_alloc(new_aM, n_new*sizeof(double), errmsg);
    class_alloc(new_dlnM2, (n_new-1)*sizeof(double), errmsg);

    for (k = 0; k < m; k++) {
      double dx = x_ext[k] - pspl->lna[0], dL = L_ext[k] - L_fit[0];
      new_x[k] = x_ext[k];
      for (j = 0; j < 3; j++)
        new_cols[j][k] = tables_smg_continued(&fit[j], pspl->y[j], _FALSE_, dx, dL);
      new_aM[k] = tables_smg_continued(&fit_aM, ptab->alpha_M[0], _FALSE_, dx, dL);
    }
    for (i = 0; i < n; i++) {
      new_x[m+i] = pspl->lna[i];
      for (j = 0; j < 3; j++) new_cols[j][m+i] = pspl->y[i*3+j];
      new_lnM2[m+i] = ptab->lnM2[i];
      new_aM[m+i] = ptab->alpha_M[i];
    }
    for (i = 0; i < n-1; i++)
      new_dlnM2[m+i] = ptab->dlnM2[i];

    /* ln M_*^2: exact integral of the continued alpha_M, backwards from the first node */
    for (k = m-1; k >= 0; k--) {
      double xa = new_x[k], xb = new_x[k+1], integral = 0.;
      for (q = 0; q < 3; q++) {
        double xq = 0.5*(xa + xb) + 0.5*(xb - xa)*gx[q], Lq, dLq;
        class_call(tables_smg_ln_rho_bg(pba, xq, &Lq, &dLq, errmsg), errmsg, errmsg);
        integral += gw[q]*tables_smg_continued(&fit_aM, ptab->alpha_M[0], _FALSE_,
                                               xq - pspl->lna[0], Lq - L_fit[0]);
      }
      new_dlnM2[k] = 0.5*(xb - xa)*integral;
      new_lnM2[k] = new_lnM2[k+1] - new_dlnM2[k];
    }

    class_call(tables_smg_free_table(pspl), errmsg, errmsg);
    free(ptab->lnM2);
    free(ptab->alpha_M);
    free(ptab->dlnM2);
    class_call(tables_smg_init_table(pspl, n_new, 3, new_x, new_cols, is_log, _FALSE_, errmsg),
               errmsg, errmsg);
    ptab->lnM2 = new_lnM2;
    ptab->alpha_M = new_aM;
    ptab->dlnM2 = new_dlnM2;
    ptab->alpha_M_slope = tables_smg_tail_slope(new_x[0], new_aM[0], new_x[1], new_aM[1], _FALSE_, _FALSE_);

    free(new_x);
    for (j = 0; j < 3; j++) free(new_cols[j]);
  }

  free(L_fit); free(x_ext); free(L_ext); free(dL_ext);
  ptab->prepared = _TRUE_;
  return _SUCCESS_;
}

/** Cubic Hermite interpolation of ln M_*^2 with slopes alpha_M. */
static void tables_smg_hermite(struct alphas_table_smg * ptab,
                               double lna,
                               double * lnM2,
                               double * dlnM2) {
  double * x = ptab->spl.lna;
  int lo = 0, hi = ptab->spl.n - 1, mid;
  double h, t, t2, t3, p0, D, m0, m1;

  while (hi - lo > 1) {
    mid = (lo + hi)/2;
    if (x[mid] > lna) hi = mid; else lo = mid;
  }
  h = x[hi] - x[lo];
  t = (lna - x[lo])/h;
  t2 = t*t;
  t3 = t2*t;
  /* standard cubic Hermite form, written with the stored increment
     D = p1 - p0 instead of p1 and p0 separately */
  p0 = ptab->lnM2[lo];
  D = ptab->dlnM2[lo];
  m0 = ptab->alpha_M[lo]*h;
  m1 = ptab->alpha_M[hi]*h;

  *lnM2 = p0 + (-2.*t3 + 3.*t2)*D + (t3 - 2.*t2 + t)*m0 + (t3 - t2)*m1;
  *dlnM2 = ((-6.*t2 + 6.*t)*D + (3.*t2 - 4.*t + 1.)*m0 + (3.*t2 - 2.*t)*m1)/h;
}

/** ln M_*^2 and alpha_M at ln a (table prepared), as in tables_smg_eval_alphas. */
static void tables_smg_lnM2(struct alphas_table_smg * ptab, double lna, double * lnM2, double * alpha_M) {
  double dlnM2;
  if (lna < ptab->spl.lna[0]) {
    double dx = lna - ptab->spl.lna[0], s = ptab->alpha_M_slope, aM0 = ptab->alpha_M[0];
    *alpha_M = aM0*exp(s*dx);
    *lnM2 = ptab->lnM2[0] + ((s != 0.) ? aM0*expm1(s*dx)/s : aM0*dx);
  }
  else {
    lna = MIN(lna, ptab->spl.lna[ptab->spl.n-1]);
    tables_smg_hermite(ptab, lna, lnM2, &dlnM2);
    *alpha_M = dlnM2;
  }
}

int tables_smg_eval_alphas(struct background * pba,
                           struct alphas_table_smg * ptab,
                           double lna,
                           double * alpha_K,
                           double * alpha_B,
                           double * alpha_T,
                           double * M2,
                           double * alpha_M,
                           ErrorMsg errmsg) {
  double values[3], lnM2;

  if (ptab->prepared == _FALSE_)
    class_call(tables_smg_prepare_alphas(pba, ptab, lna, errmsg), errmsg, errmsg);

  class_call(tables_smg_eval_table(&(ptab->spl), lna, values, NULL, errmsg), errmsg, errmsg);
  *alpha_K = values[0];
  *alpha_B = values[1];
  *alpha_T = values[2];

  /* below the first node: alpha_M = alpha_M0 exp(s dx) (or constant), and
     ln M_*^2 = ln M_*^2_0 + int_0^dx alpha_M, exactly */
  tables_smg_lnM2(ptab, lna, &lnM2, alpha_M);
  *M2 = exp(lnM2);
  return _SUCCESS_;
}

int tables_smg_free_alphas(struct alphas_table_smg * ptab) {
  tables_smg_free_table(&(ptab->spl));
  free(ptab->lnM2);
  free(ptab->alpha_M);
  free(ptab->dlnM2);
  return _SUCCESS_;
}

/* ------------------------------------------------------------------------ */
/* expansion table                                                          */
/* ------------------------------------------------------------------------ */

int tables_smg_read_expansion(struct file_content * pfc,
                              struct expansion_table_smg * ptab,
                              ErrorMsg errmsg) {

  /* column order: ln a, rho_smg/rho_smg(a=1) [, w_smg] */
  char * keys[3] = {"table_lna_de_smg", "table_rho_de_smg", "table_w_de_smg"};
  double * cols[3];
  double * spl_cols[2];
  short is_log[2] = {_TRUE_, _FALSE_};
  short constant;
  int n, n_found, i, j, flag;
  double tol = 1e-6, lnrho0;

  class_call(tables_smg_read_columns(pfc, "expansion_table_file_smg", keys, 3, 2, &n, &n_found, cols, errmsg),
             errmsg, errmsg);
  class_call(tables_smg_validate("expansion table", n, n_found, cols, errmsg), errmsg, errmsg);
  class_call(tables_smg_read_extrapolation(pfc, &(ptab->extrapolation), errmsg), errmsg, errmsg);
  class_call(parser_read_double(pfc, "table_rho_norm_tol_smg", &tol, &flag, errmsg), errmsg, errmsg);
  if (flag == _FALSE_) tol = 1e-6;

  for (i = 0; i < n; i++) {
    class_test(cols[1][i] <= 0., errmsg,
               "expansion table: rho_smg must be positive, got %g at node %d", cols[1][i], i);
    cols[1][i] = log(cols[1][i]);
  }

  constant = (ptab->extrapolation.mode == table_constant) ? _TRUE_ : _FALSE_;
  ptab->has_w = (n_found == 3) ? _TRUE_ : _FALSE_;
  spl_cols[0] = cols[1];
  spl_cols[1] = cols[2];
  class_call(tables_smg_init_table(&(ptab->spl), n, n_found - 1, cols[0], spl_cols, is_log, constant, errmsg),
             errmsg, errmsg);

  /* normalise exactly at a = 1, so that rho_smg(a=1) = Omega0_smg H0^2 and H(a=1) = H0;
     a large mismatch means the table is not normalised (e.g. given in physical units) */
  {
    double values[2];
    class_call(tables_smg_eval_table(&(ptab->spl), 0., values, NULL, errmsg), errmsg, errmsg);
    lnrho0 = values[0];
  }
  class_test(fabs(expm1(lnrho0)) > tol, errmsg,
             "expansion table: rho_smg(a=1) = %.10g, but the table must be normalised to 1 at a = 1 "
             "(tolerance table_rho_norm_tol_smg = %g)", exp(lnrho0), tol);
  for (i = 0; i < n; i++)
    ptab->spl.y[i*ptab->spl.ncol] -= lnrho0; /* spline of ln rho is linear in the data */

  /* the scaling/power-law continuation needs hi_class's background: done at first use */
  ptab->prepared = constant;

  for (j = 0; j < 3; j++) free(cols[j]);
  return _SUCCESS_;
}

/** Add the early-time continuation of the expansion history (see tables_smg.h). */
static int tables_smg_prepare_expansion(struct background * pba,
                                        struct expansion_table_smg * ptab,
                                        double lna_check,
                                        ErrorMsg errmsg) {
  struct table_smg * pspl = &(ptab->spl);
  int n = pspl->n, ncol = pspl->ncol, n_fit, m, k, i, j, n_new;
  double * L_fit, * x_ext, * L_ext, * dL_ext, * new_x, * new_cols[2] = {NULL, NULL};
  struct table_scaling_fit_smg fit;
  short is_log[2] = {_TRUE_, _FALSE_};

  class_call(tables_smg_continuation_grid(pba, &(ptab->extrapolation), n, pspl->lna,
                                          &n_fit, &L_fit, &m, &x_ext, &L_ext, &dL_ext, errmsg),
             errmsg, errmsg);

  if (m > 0) {
    double L_check, dummy, rho_unit = pba->Omega0_smg*pba->H0*pba->H0;
    short use_field = _FALSE_;
    double * lnF = NULL, F0 = 0.;

    class_call(tables_smg_ln_rho_bg(pba, lna_check, &L_check, &dummy, errmsg), errmsg, errmsg);

    /* With tabulated alphas, continue the field's own energy density instead of
       rho_smg: 3 M_*^2 H^2 = rho_bg + rho_phi, i.e.
         rho_phi = M_*^2 (rho_bg + rho_smg) - rho_bg,
       and rebuild rho_smg = (rho_bg + rho_phi)/M_*^2 - rho_bg with the continued
       M_*^2 of the alphas table. This keeps the expansion history consistent with
       M_*^2 (e.g. Brans-Dicke, where rho_smg ~ (1/M_*^2 - 1) rho_bg early on). */
    if ((pba->gravity_model_smg == alphas_table) && (pba->alphas_table_smg != NULL) && (rho_unit > 0.)) {
      double lnM2, aM, Fi;
      short ok = _TRUE_;
      /* prepare the alphas here, so that their reliability check uses the earliest a */
      if (pba->alphas_table_smg->prepared == _FALSE_)
        class_call(tables_smg_prepare_alphas(pba, pba->alphas_table_smg, lna_check, errmsg), errmsg, errmsg);
      class_alloc(lnF, n_fit*sizeof(double), errmsg);
      for (i = 0; i < n_fit; i++) {
        /* rho_phi = (M_*^2 - 1) rho_bg + M_*^2 rho_smg, free of cancellations */
        tables_smg_lnM2(pba->alphas_table_smg, pspl->lna[i], &lnM2, &aM);
        Fi = expm1(lnM2)*exp(L_fit[i]) + exp(lnM2)*rho_unit*exp(pspl->y[i*ncol]);
        if (i == 0) F0 = Fi;
        if ((Fi == 0.) || ((Fi > 0.) != (F0 > 0.))) ok = _FALSE_;
        lnF[i] = (Fi != 0.) ? log(fabs(Fi)) : 0.;
      }
      use_field = ok;
    }

    /* direct continuation of rho_smg */
    fit = tables_smg_fit_scaling(&(ptab->extrapolation), n_fit, pspl->lna, L_fit, pspl->y, ncol, _TRUE_);

    /* the field route is used when the scaling family describes rho_phi better
       than rho_smg, and when it keeps rho_smg > 0 at all continuation nodes */
    if (use_field == _TRUE_) {
      struct table_scaling_fit_smg fit_F;
      fit_F = tables_smg_fit_scaling(&(ptab->extrapolation), n_fit, pspl->lna, L_fit, lnF, 1, _TRUE_);
      if ((fit_F.constant == _FALSE_) && (fit_F.rms < fit.rms)) {
        for (k = 0; k < m; k++) {
          double lnM2, aM, F;
          tables_smg_lnM2(pba->alphas_table_smg, x_ext[k], &lnM2, &aM);
          F = F0*exp(fit_F.s*(x_ext[k] - pspl->lna[0]) + fit_F.beta*(L_ext[k] - L_fit[0]));
          if (F - expm1(lnM2)*exp(L_ext[k]) <= 0.) { use_field = _FALSE_; break; }
        }
      }
      else {
        use_field = _FALSE_;
      }
      if (use_field == _TRUE_) fit = fit_F;
      if (ptab->extrapolation.verbose > 1)
        printf(" -> expansion table continued through %s (fit rms: rho_smg %.2e, rho_phi %.2e)\n",
               use_field == _TRUE_ ? "rho_phi = M_*^2 (rho_bg + rho_smg) - rho_bg" : "rho_smg",
               tables_smg_fit_scaling(&(ptab->extrapolation), n_fit, pspl->lna, L_fit, pspl->y, ncol, _TRUE_).rms,
               fit_F.rms);
    }

    if (use_field == _TRUE_) {
      class_call(tables_smg_check_fit(&(ptab->extrapolation), "rho_phi", &fit, n_fit, pspl->lna, L_fit,
                                      lnF, 1, _TRUE_, lna_check, L_check, errmsg),
                 errmsg, errmsg);
    }
    else {
      class_call(tables_smg_check_fit(&(ptab->extrapolation), "rho_smg", &fit, n, pspl->lna, L_fit,
                                      pspl->y, ncol, _TRUE_, lna_check, L_check, errmsg),
                 errmsg, errmsg);
    }

    n_new = m + n;
    class_alloc(new_x, n_new*sizeof(double), errmsg);
    for (j = 0; j < ncol; j++) class_alloc(new_cols[j], n_new*sizeof(double), errmsg);

    for (k = 0; k < m; k++) {
      double dx = x_ext[k] - pspl->lna[0], dL = L_ext[k] - L_fit[0];
      new_x[k] = x_ext[k];
      if (use_field == _TRUE_) {
        double lnM2, aM, M2, D, rho_bg = exp(L_ext[k]), F, dF, rho_s, drho_s;
        tables_smg_lnM2(pba->alphas_table_smg, x_ext[k], &lnM2, &aM);
        M2 = exp(lnM2);
        D = expm1(lnM2);
        F = (fit.constant == _TRUE_) ? F0 : F0*exp(fit.s*dx + fit.beta*dL);
        dF = (fit.constant == _TRUE_) ? 0. : F*(fit.s + fit.beta*dL_ext[k]);
        /* rho_smg = (rho_phi - (M_*^2 - 1) rho_bg)/M_*^2, free of cancellations */
        rho_s = (F - D*rho_bg)/M2;
        class_test(rho_s <= 0., errmsg,
                   "expansion table: the early-time continuation gives rho_smg <= 0 at a = %e "
                   "(rho_phi = %e, M_*^2 = %e); start the table earlier", exp(x_ext[k]), F, M2);
        /* d rho_smg / d ln a, with d ln M_*^2 / d ln a = alpha_M */
        drho_s = (dF - M2*aM*rho_bg - D*rho_bg*dL_ext[k])/M2 - aM*rho_s;
        new_cols[0][k] = log(rho_s/rho_unit);
        if (ncol > 1) new_cols[1][k] = -1. - drho_s/(3.*rho_s);
      }
      else {
        new_cols[0][k] = tables_smg_continued(&fit, pspl->y[0], _TRUE_, dx, dL);
        if (ncol > 1)
          new_cols[1][k] = (fit.constant == _TRUE_) ? -1. : -1. - (fit.s + fit.beta*dL_ext[k])/3.;
      }
    }
    free(lnF);
    for (i = 0; i < n; i++) {
      new_x[m+i] = pspl->lna[i];
      for (j = 0; j < ncol; j++) new_cols[j][m+i] = pspl->y[i*ncol+j];
    }

    class_call(tables_smg_free_table(pspl), errmsg, errmsg);
    class_call(tables_smg_init_table(pspl, n_new, ncol, new_x, new_cols, is_log, _FALSE_, errmsg),
               errmsg, errmsg);

    free(new_x);
    for (j = 0; j < ncol; j++) free(new_cols[j]);
  }

  free(L_fit); free(x_ext); free(L_ext); free(dL_ext);
  ptab->prepared = _TRUE_;
  return _SUCCESS_;
}

int tables_smg_eval_expansion(struct background * pba,
                              struct expansion_table_smg * ptab,
                              double lna,
                              double * rho_over_rho0,
                              double * w,
                              ErrorMsg errmsg) {
  double values[2], dvalues[2];

  if (ptab->prepared == _FALSE_)
    class_call(tables_smg_prepare_expansion(pba, ptab, lna, errmsg), errmsg, errmsg);

  class_call(tables_smg_eval_table(&(ptab->spl), lna, values, dvalues, errmsg), errmsg, errmsg);
  *rho_over_rho0 = exp(values[0]);

  if ((ptab->has_w == _TRUE_) && (lna >= ptab->spl.lna[0])) {
    *w = values[1];
  }
  else {
    /* from the continuity equation, d ln rho/d ln a = -3 (1 + w); below the
       first node this keeps w consistent with the extrapolated rho */
    *w = -1. - dvalues[0]/3.;
  }
  return _SUCCESS_;
}

int tables_smg_free_expansion(struct expansion_table_smg * ptab) {
  tables_smg_free_table(&(ptab->spl));
  return _SUCCESS_;
}
