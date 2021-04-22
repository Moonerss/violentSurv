#' @title Stepwise Variable Selection Procedure for Cox's Proportional Hazards Model and Cox's Model
#'
#' @description
#' This stepwise variable selection procedure (with iterations between the 'forward' and 'backward' steps) can be applied to obtain the best candidate final Cox's proportional hazards model or Cox's proportional hazards model with time-dependent covariates (called the Cox's model).
#' @details The goal of regression analysis is to find one or a few parsimonious regression models that fit the observed data well for effect estimation and/or outcome prediction. To ensure a good quality of analysis, the model-fitting techniques for (1) variable selection, (2) goodness-of-fit assessment, and (3) regression diagnostics and remedies should be used in regression analysis. The stepwise variable selection procedure (with iterations between the 'forward' and 'backward' steps) is one of the best ways to obtaining the best candidate final regression model. All the bivariate significant and non-significant relevant covariates and some of their interaction terms (or moderators) are put on the 'variable list' to be selected. The significance levels for entry (SLE) and for stay (SLS) are suggested to be set at 0.15 or larger for being conservative. Then, with the aid of substantive knowledge, the best candidate final regression model is identified manually by dropping the covariates with p value > 0.05 one at a time until all regression coefficients are significantly different from 0 at the chosen alpha level of 0.05. Since the statistical testing at each step of the stepwise variable selection procedure is conditioning on the other covariates in the regression model, the multiple testing problem is not of concern. Any discrepancy between the results of bivariate analysis and regression analysis is likely due to the confounding effects of uncontrolled covariates in bivariate analysis or the masking effects of intermediate variables (or mediators) in regression analysis.
#' @param Time The 'Time' (time to an event) for the sepcified Cox's proportional hazards model as in coxph().
#' @param T1 The 'T1' (Start) of the long-form data for the sepcified Cox's model as in coxph().
#' @param T2 The 'T2' (Stop) of the long-form data for the sepcified Cox's model as in coxph().
#' @param Status The 'Status' (event indicator) for the sepcified Cox's proportional hazards model as in coxph().
#' @param variable.list A list of covariates to be selected.
#' @param in.variable A list of covariate(s) to be always included in the regression model.
#' @param data The data to be analyzed.
#' @param sle The chosen significance level for entry (SLE).
#' @param sls The chosen significance level for stay (SLS).
#' @param vif.threshold The chosen threshold value of variance inflating factor (VIF).
#' @return
#' A model object representing the identified "Stepwise Final Model" with the values of variance inflating factor (VIF) for all included covarites is displayed.
#' @examples
#' \dontrun{
#' The data 'lung' is available in the 'survival' package.}
#'
#'if (requireNamespace("survival", quietly = TRUE)) {
#'  lung <- survival::lung
#'}
#'
#'names(lung)
#'dim(lung)
#'my.data <- na.omit(lung)
#'dim(my.data)
#'head(my.data)
#'my.data$status1 <- ifelse(my.data$status==2,1,0)
#'my.variable.list <- c("inst", "age", "sex", "ph.ecog", "ph.karno", "pat.karno")
#'My.stepwise.coxph(Time = "time", Status = "status1", variable.list = my.variable.list,
#'    in.variable = c("meal.cal", "wt.loss"), data = my.data)
#'
#'my.variable.list <- c("inst", "age", "sex", "ph.ecog", "ph.karno", "pat.karno", "meal.cal",
#'    "wt.loss")
#'My.stepwise.coxph(Time = "time", Status = "status1", variable.list = my.variable.list,
#'    data = my.data, sle = 0.25, sls = 0.25)
#' @section Warning:
#' The value of variance inflating factor (VIF) is bigger than 10 in continuous covariates or VIF is bigger than 2.5 in categorical covariates indicate the occurrence of multicollinearity problem among some of the covariates in the fitted regression model.
#' @export
#' @import car
#' @import survival
#' @importFrom stats anova binomial glm update

stepwise_coxph <- function(Time=NULL, T1=NULL, T2=NULL, Status=NULL, variable.list, in.variable="NULL", data, sle=0.15, sls=0.15, vif.threshold=999)
{

  univar.pvalue <- NULL
  temp.model <- NULL

  if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time,", ", Status,") ~ ", paste(in.variable, collapse="+"), sep="")), data=data, method="efron")
  } else if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time,", ", Status,") ~ ", paste(in.variable, collapse="+"), sep="")), data=data, method="efron")
  } else if (is.null(Time)){
    initial.model <- coxph(as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", paste(in.variable, collapse="+"), sep="")), data=data, method="efron")
  } else if (is.null(Time)){
    initial.model <- coxph(as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", paste(in.variable, collapse="+"), sep="")), data=data, method="efron")
  }

  if (is.null(initial.model$coefficients)) {

    for (i in 1:length(variable.list))
    {
      if (is.null(T2))
      {
        uni.model <- coxph(as.formula(paste("Surv(", Time,", ", Status,") ~ ", variable.list[i], sep="")), data=data, method="efron")
      }
      if (is.null(Time))
      {
        uni.model <- coxph(as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", variable.list[i], sep="")), data=data, method="efron")
      }
      univar.pvalue[i] <- summary(uni.model)$coefficients[5]
    }

    variable.list1 <- variable.list[univar.pvalue<=0.9 & !is.na(univar.pvalue)]
    univar.pvalue1 <- univar.pvalue[univar.pvalue<=0.9 & !is.na(univar.pvalue)]
    uni.x <- variable.list1[which.min(univar.pvalue1)]
    if (length(uni.x) > 0) {
      if (is.null(T2))
      {
        formula <- as.formula(paste("Surv(", Time,", ", Status,") ~ ", uni.x, sep=""))
        temp.model <- coxph(formula, data=data, method="efron")
        if (length(temp.model$coefficients) > 1)
        {
          print(vif(glm(paste(Status, paste(names(temp.model$coefficients), collapse="+"), sep="~"), data=data, family=binomial(link="logit"))))
        }
      }
      if (is.null(Time))
      {
        formula <- as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", uni.x, sep=""))
        temp.model <- coxph(formula, data=data, method="efron")
      }

      cat("# --------------------------------------------------------------------------------------------------\n")
      cat("# Initial Model:\n")
      print(summary(temp.model))

    }

  } else if (!is.null(initial.model$coefficients)) {
    temp.model <- initial.model
    cat("# --------------------------------------------------------------------------------------------------\n")
    cat("# Initial Model:\n")
    print(summary(temp.model))

  }
  if (length(temp.model$coefficients) > 1)
  {
    cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
    cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
    print(vif(glm(paste(Status, paste(names(temp.model$coefficients), collapse="+"), sep="~"), data=data, family=binomial(link="logit"))))
  }

  i <- 0
  break.rule <- TRUE
  while (break.rule)
  {
    i <- i + 1
    if (i == 1)
    {
      variable.list2 <- setdiff(variable.list, all.vars(temp.model$formula))
    } else
    {
      variable.list2 <- setdiff(variable.list, c(all.vars(temp.model$formula), out.x))
      out.x <- NULL
    }

    if (length(variable.list2) != 0)
    {
      anova.pvalue <- NULL
      mv.pvalue <- NULL
      vif.value <- NULL
      for (k in 1:length(variable.list2))
      {
        model <- update(temp.model, as.formula(paste(". ~ . + ", variable.list2[k], sep="")))
        if (length(model$coefficients) > 1)
        {
          if (sum(is.na(model$coefficients)) != 0)
          {
            anova.pvalue[k] <- 1
            mv.pvalue[k] <- 1
            vif.value[k] <- 999
          } else {
            anova.pvalue[k] <- anova(temp.model, model)[2,"P(>|Chi|)"]
            mv.pvalue[k] <- summary(model)$coefficients[nrow(summary(model)$coefficients),"Pr(>|z|)"]
            model.vif <- vif(glm(as.formula(paste(Status, paste(names(model$coefficients), collapse="+"), sep="~")), data=data, family=binomial(link="logit")))
            vif.value[k] <- model.vif[length(model.vif)]
          }
        }
      }

      variable.list2.1 <- variable.list2[mv.pvalue<=0.9 & !is.na(mv.pvalue) & vif.value <= vif.threshold]
      anova.pvalue2 <- anova.pvalue[mv.pvalue<=0.9 & !is.na(mv.pvalue) & vif.value <= vif.threshold]
      mv.pvalue2 <- mv.pvalue[mv.pvalue<=0.9 & !is.na(mv.pvalue) & vif.value <= vif.threshold]
      enter.x <- variable.list2.1[anova.pvalue2==min(anova.pvalue2, na.rm=TRUE) & anova.pvalue2 <= sle]
      wald.p <- mv.pvalue2[anova.pvalue2==min(anova.pvalue2, na.rm=TRUE) & anova.pvalue2 <= sle]
      if (length(setdiff(enter.x, NA)) != 0)
      {
        if (length(enter.x) > 1)
        {
          enter.x <- enter.x[which.min(wald.p)]
        }

        cat("# --------------------------------------------------------------------------------------------------", "\n")
        cat(paste("### iter num = ", i, ", Forward Selection by LR Test: ","+ ", enter.x, sep=""), "\n")
        temp.model <- update(temp.model, as.formula(paste(". ~ . + ", enter.x, sep="")))
        print(summary(temp.model))
        cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
        cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
        if (length(temp.model$coefficients) > 1)
        {
          print(vif(glm(paste(Status, paste(names(temp.model$coefficients), collapse="+"), sep="~"), data=data, family=binomial(link="logit"))))
        }
      }
    } else {enter.x <- NULL}

    if (i == 1 & length(enter.x) == 0) {
      cat("# ==================================================================================================", "\n")
      cat(paste("*** Stepwise Final Model (in.lr.test: sle = ", sle, "; variable selection restrict in vif = ", vif.threshold, "):", sep=""), "\n")
      print(summary(temp.model))
      break
    } else {
      variable.list3 <- setdiff(rownames(summary(temp.model)$coefficients), c(enter.x, in.variable))
      if (length(variable.list3) != 0)
      {
        anova.pvalue <- NULL
        for (k in 1:length(variable.list3))
        {
          model <- update(temp.model, as.formula(paste(". ~ . - ", variable.list3[k], sep="")))
          anova.pvalue[k] <- anova(model, temp.model)[2,"P(>|Chi|)"]
        }

        out.x <- variable.list3[anova.pvalue==max(anova.pvalue, na.rm=TRUE) & anova.pvalue > sls]
        out.x <- setdiff(out.x, NA)
        if (length(out.x) != 0)
        {
          if (length(out.x) > 1)
          {
            out.x.1 <- out.x
            for (j in 1:length(out.x)) {
              out.x[j] <- out.x.1[(length(out.x)-j+1)]
            }

            wald.p <- rep(NA, length(out.x))
            for (j in 1:length(out.x)) {
              wald.p[j] <- summary(temp.model)$coefficients[,"Pr(>|z|)"][rownames(summary(temp.model)$coefficients)==out.x[j]]
            }
            out.x <- out.x[which.max(wald.p)]
          }
          cat("# --------------------------------------------------------------------------------------------------", "\n")
          cat(paste("### iter num = ", i, ", Backward Selection by LR Test: ","- ", out.x, sep=""), "\n")
          temp.model <- update(temp.model, as.formula(paste(". ~ . - ", out.x, sep="")))
          print(summary(temp.model))
          cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
          cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
          if (length(temp.model$coefficients) > 1)
          {
            print(vif(glm(paste(Status, paste(names(temp.model$coefficients), collapse="+"), sep="~"), data=data, family=binomial(link="logit"))))
          }
        }
      } else {out.x <- NULL}
    }

    if ((length(enter.x) + length(out.x)) == 0)
    {
      final.model <- temp.model
      cat("# ==================================================================================================", "\n")
      cat(paste("*** Stepwise Final Model (in.lr.test: sle = ", sle, "; out.lr.test: sls = ", sls, "; variable selection restrict in vif = ", vif.threshold, "):", sep=""), "\n")
      print(summary(final.model))
      cat("--------------- Variance Inflating Factor (VIF) ---------------\n")
      cat("Multicollinearity Problem: Variance Inflating Factor (VIF) is bigger than 10 (Continuous Variable) or is bigger than 2.5 (Categorical Variable)\n")
      if (length(final.model$coefficients) > 1)
      {
        print(vif(glm(paste(Status, paste(names(final.model$coefficients), collapse="+"), sep="~"), data=data, family=binomial(link="logit"))))
      }
      break.rule <- FALSE
    }

    enter.x <- NULL
  }
}
