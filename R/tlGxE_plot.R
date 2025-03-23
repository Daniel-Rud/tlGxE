

#'
#' Manhattan Plots for ATE and MOR additive and codominant tests.
#'
#' @param object a \code{tlGxE} object, ie, a reference to the \code{tlGxE} function call.
#' @param ... see details
#' @description
#' Manhattan plots for ATE and MOR additive and codominant tests in \code{tlGxE}.
#' @details
#' Pass the desired testing significance level as an argument in the function call.  Default level is \eqn{\alpha = 0.05}, which should
#' be changed to account for multiple testing.
#' @export

plot.tlGxE = function(object, ...)
{
  # check for fgwqsr object
  if(!methods::is(object, "tlGxE"))
  {
    stop("Must pass the object output of the tlGxE() function.")
  }

  passed_args = list(...)

  alpha = 0.05
  if(length(passed_args) >= 1)
  {
    alpha = passed_args[[1]]
  }

  sig_lvl = -1*log(alpha, base = 10)

  result_frame = object
  G_names = colnames(result_frame)

  if(length(passed_args) == 1)
  {
    sig_lvl = -1*log(passed_args[[1]], base = 10)
  }

  #############################################################################
  # ATE Plots #################################################################
  #############################################################################

  # ATE estimate plot


  # ATE pvalue plot

  ATE_aov_pvalues = result_frame["ATE_codominant_pvalue", ]
  ATE_lin_pvalues = result_frame["ATE_additive_pvalue", ]

  ATE_signif_aov = ifelse(-1*log(ATE_aov_pvalues, base = 10) < sig_lvl, 0, 1) %>% factor
  ATE_signif_lin = ifelse(-1*log(ATE_lin_pvalues, base = 10) < sig_lvl, 0, 1) %>% factor

  ATE_aov_data = data.frame(G_names = G_names %>% factor(levels =G_names, ordered = T),
                            neg_log_10_pvalues = -1 * log(ATE_aov_pvalues, base = 10),
                            Significant = ATE_signif_aov)
  ATE_linear_data = data.frame(G_names = G_names %>% factor(levels =G_names, ordered = T),
                               neg_log_10_pvalues = -1 * log(ATE_lin_pvalues, base = 10),
                               Significant = ATE_signif_lin)

  ATE_aov_plot = ggplot2::ggplot(data = ATE_aov_data, mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = Significant)) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::geom_hline(yintercept = sig_lvl, color = "red", linetype = "dashed") +
    ggplot2::ylab(latex2exp::TeX("$-log_{10}(pvalue)$")) + ggplot2::xlab("SNP") +
    ggplot2::geom_text(mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues, label =G_names),
                       data = subset(ATE_aov_data, ATE_aov_data$Significant == 1),
                       check_overlap = F, hjust = -0.2, vjust = .5, size = 3.5) +
    ggplot2::ggtitle("tlGxE ATE Codominant") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 2))

  ATE_linear_plot = ggplot2::ggplot(data = ATE_linear_data, mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = Significant)) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::geom_hline(yintercept = sig_lvl, color = "red", linetype = "dashed") +
    ggplot2::ylab(latex2exp::TeX("$-log_{10}(pvalue)$")) + ggplot2::xlab("SNP") +
    ggplot2::geom_text(mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues, label =G_names),
                       data = subset(ATE_linear_data, ATE_linear_data$Significant == 1),
                       check_overlap = F, hjust = -0.2, vjust = .5, size = 3.5) +
    ggplot2::ggtitle("tlGxE ATE Additive") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 2))




  #############################################################################

  # check if MOR plots should be included
  plot_mors = "MOR_G0" %in% rownames(result_frame)

  if(plot_mors)
  {
    #############################################################################
    # MOR Plots #################################################################
    #############################################################################

    # MOR estimate plot
    # MOR pvalue plot

    MOR_aov_pvalues = result_frame["MOR_codominant_pvalue", ]
    MOR_mult_pvalues = result_frame["MOR_additive_pvalue", ]

    MOR_signif_aov = ifelse(-1*log(MOR_aov_pvalues, base = 10) < sig_lvl, 0, 1) %>% factor
    MOR_signif_mult = ifelse(-1*log(MOR_mult_pvalues, base = 10) < sig_lvl, 0, 1) %>% factor



    MOR_aov_data = data.frame(G_names = G_names %>% factor(levels =G_names, ordered = T),
                              neg_log_10_pvalues = -1*log(MOR_aov_pvalues, base = 10),
                              Significant = MOR_signif_aov)
    MOR_mult_data = data.frame(G_names = G_names %>% factor(levels =G_names, ordered = T),
                               neg_log_10_pvalues = -1*log(MOR_mult_pvalues, base = 10),
                               Significant = MOR_signif_mult)

    MOR_aov_plot = ggplot2::ggplot(data = MOR_aov_data, mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues)) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = Significant)) +
      ggplot2::scale_color_manual(values = c("black", "red")) +
      ggplot2::geom_hline(yintercept = sig_lvl, color = "red", linetype = "dashed") +
      ggplot2::ylab(latex2exp::TeX("$-log_{10}(pvalue)$")) + ggplot2::xlab("SNP") +
      ggplot2::geom_text(mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues, label =G_names),
                         data = subset(MOR_aov_data, MOR_aov_data$Significant == 1),
                         check_overlap = F,hjust = -0.2, vjust = .5, size = 3.5) +
      ggplot2::ggtitle("tlGxE MOR Codominant") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::guides(color = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 2))


    MOR_mult_plot = ggplot2::ggplot(data = MOR_mult_data, mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues)) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = Significant)) +
      ggplot2::scale_color_manual(values = c("black", "red")) +
      ggplot2::geom_hline(yintercept = sig_lvl, color = "red", linetype = "dashed") +
      ggplot2::ylab(latex2exp::TeX("$-log_{10}(pvalue)$")) + ggplot2::xlab("SNP") +
      ggplot2::geom_text(mapping = ggplot2::aes(x = G_names, y = neg_log_10_pvalues, label =G_names),
                         data = subset(MOR_mult_data, MOR_mult_data$Significant == 1),
                         check_overlap = F,hjust = -0.2, vjust = .5, size = 3.5) +
      ggplot2::ggtitle("tlGxE MOR Additive") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::guides(color = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 2))
  }

  #############################################################################
  if(plot_mors)
  {
    plot = ggpubr::ggarrange(ATE_aov_plot, ATE_linear_plot,
                             MOR_aov_plot, MOR_mult_plot,
                             nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
  }else
  {
    plot = ggpubr::ggarrange(ATE_aov_plot, ATE_linear_plot,
                             nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
  }

  return(plot)
}
