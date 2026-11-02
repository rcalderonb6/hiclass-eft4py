#ifndef __TABLES_SMG__
#define __TABLES_SMG__

#include "common.h"
#include "parser.h"
#include "arrays.h"

/**
 * Tabulated inputs for hi_class parametrisations (_smg):
 *  - gravity_model = alphas_table : alpha_K, alpha_B, alpha_M, alpha_T, M_*^2 as functions of ln a
 *  - expansion_model = table      : rho_smg(a)/rho_smg(a=1), optionally w_smg(a)
 *
 * Tables are given either as a file (columns separated by blanks, '#' comments
 * allowed) or as comma-separated lists of equal length (from the Python wrapper;
 * .ini lines are limited to _LINE_LENGTH_MAX_ characters).
 *
 * Interpolation is a cubic spline in ln a. M_*^2 and alpha_M are interpolated
 * together with a cubic Hermite interpolant of ln M_*^2 whose node slopes are
 * alpha_M, so that alpha_M = d ln M_*^2 / d ln a holds exactly everywhere.
 * rho_smg is splined in ln rho_smg.
 *
 * Before the first node (table_extrapolation_smg):
 *  - scaling (default): each of alpha_K, alpha_B, alpha_T, alpha_M and rho_smg
 *    is continued as f = f_0 (a/a_0)^s (rho_bg(a)/rho_bg(a_0))^beta, where
 *    rho_bg is hi_class's own radiation + matter density (photons, ur, baryons,
 *    cdm, ncdm). s and beta are a least-squares fit over the first
 *    table_extrapolation_window_smg e-folds of nodes, pinned to the first node.
 *    This family is exact for constants, power laws, quantities tracking the
 *    background (rho_smg propto rho_bg) and quantities sourced by the
 *    matter/radiation ratio (e.g. alpha propto rho_m/rho_bg), across
 *    matter-radiation equality. When the window is too deep in radiation
 *    domination to tell s from beta, beta = 0 (equivalent there).
 *  - power_law: the same with beta = 0.
 *  - constant: every column held at its first-node value.
 *  A column that vanishes or changes sign within the window is held constant.
 *  Reliability check: the fit is repeated over half the window; if the two
 *  continuations of any column differ by more than table_extrapolation_tol_smg
 *  (default 0.2) in ln f at the earliest a hi_class uses, the table does not
 *  determine the early-time behaviour and an error asks for an earlier start.
 *  ln M_*^2 is the exact integral of the extrapolated alpha_M, and w follows
 *  from the extrapolated rho_smg. With tabulated alphas, rho_smg may instead be
 *  continued through the field's own energy density, from the Friedmann
 *  equation 3 M_*^2 H^2 = rho_bg + rho_phi: rho_phi = (M_*^2 - 1) rho_bg +
 *  M_*^2 rho_smg is continued and rho_smg rebuilt with the continued M_*^2, so
 *  that the expansion stays consistent with M_*^2 (e.g. Brans-Dicke, where
 *  rho_smg ~ (1/M_*^2 - 1) rho_bg). This route is taken when the scaling family
 *  fits rho_phi better than rho_smg over the window (and keeps rho_smg > 0). The continuation is tabulated (prepended to
 *  the table) down to a = 1e-30 the first time the table is evaluated, and
 *  continued below that with the local power law.
 *
 * Evaluation after the last node is an error: tables must reach a >= 1.
 */

struct background;

enum table_extrapolation_smg {table_scaling, table_power_law, table_constant};

/** One table: n nodes in ln a, ncol columns (line-major, as in CLASS arrays). */
struct table_smg {
  int n;                 /**< number of nodes */
  int ncol;              /**< number of columns, excluding ln a */
  double * lna;          /**< ln a at the nodes, strictly increasing */
  double * y;            /**< values, y[i*ncol+j] */
  double * ddy;          /**< spline second derivatives, same layout */
  short * is_log;        /**< column j is a logarithm: extrapolated linearly instead of as a power law */
  double * extrap_slope; /**< below the first node: d y_j/d ln a if is_log[j], else d ln|y_j|/d ln a (0 = constant) */
};

/** Extrapolation settings, shared by both tables. */
struct table_extrapolation_settings_smg {
  enum table_extrapolation_smg mode;
  double window;         /**< fit window, in e-folds after the first node */
  double tol;            /**< reliability check: maximum |ln f_window - ln f_half window| at the
                              earliest a used (0 or negative: no check) */
  int verbose;           /**< background_verbose, for diagnostics of the continuation */
};

/** Tabulated alphas: spline columns and the Hermite data for M_*^2. */
struct alphas_table_smg {
  struct table_smg spl;  /**< columns: alpha_K, alpha_B, alpha_T */
  double * lnM2;         /**< ln M_*^2 at the nodes */
  double * alpha_M;      /**< alpha_M at the nodes = Hermite slopes of ln M_*^2 */
  double * dlnM2;        /**< ln M_*^2 increment over each interval (n-1 values), stored rather than
                              differenced so that tiny alpha_M on top of M_*^2 != 1 do not cancel out */
  double alpha_M_slope;  /**< d ln|alpha_M| / d ln a below the first node (0 = constant);
                              M_*^2 follows by integrating alpha_M exactly */
  struct table_extrapolation_settings_smg extrapolation;
  short prepared;        /**< whether the early-time continuation has been added */
};

/** Tabulated expansion history. */
struct expansion_table_smg {
  struct table_smg spl;  /**< columns: ln(rho/rho_0) [, w] */
  short has_w;           /**< whether w is tabulated (then p = w rho exactly) */
  struct table_extrapolation_settings_smg extrapolation;
  short prepared;        /**< whether the early-time continuation has been added */
};

#ifdef __cplusplus
extern "C" {
#endif

  int tables_smg_read_alphas(struct file_content * pfc,
                             struct alphas_table_smg * ptab,
                             ErrorMsg errmsg);

  int tables_smg_read_expansion(struct file_content * pfc,
                                struct expansion_table_smg * ptab,
                                ErrorMsg errmsg);

  int tables_smg_eval_alphas(struct background * pba,
                             struct alphas_table_smg * ptab,
                             double lna,
                             double * alpha_K,
                             double * alpha_B,
                             double * alpha_T,
                             double * M2,
                             double * alpha_M,
                             ErrorMsg errmsg);

  int tables_smg_eval_expansion(struct background * pba,
                                struct expansion_table_smg * ptab,
                                double lna,
                                double * rho_over_rho0,
                                double * w,
                                ErrorMsg errmsg);

  int tables_smg_free_alphas(struct alphas_table_smg * ptab);

  int tables_smg_free_expansion(struct expansion_table_smg * ptab);

#ifdef __cplusplus
}
#endif

#endif
