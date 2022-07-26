## ----setup1, include=FALSE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.asp=0.625, fig.path = 'figs/', message = FALSE,  warnings = FALSE,  cache = TRUE, eval=TRUE, fig.align='center', fig.show='hold')
options(digits=4)


## ----reading------------------------------------------------------------------------------------------------------------
library(nlme)
library(MASS)
library(lattice)
pisa <- read.table(file = "data/pisaSEL.txt", header = TRUE) 
pisa$CESCS <- pisa$ESCS - pisa$ESCSM
pisa$SEX <- factor(pisa$SEX)
pisa$CHANGE <- factor(pisa$CHANGE)
summary(pisa)


## ----grouped------------------------------------------------------------------------------------------------------------
pisa <- groupedData(MATH ~ ESCS | SCHOOLID, data = pisa)
table(getGroups(pisa))    


## ----scatter1-----------------------------------------------------------------------------------------------------------
plot(MATH ~ ESCS, pisa, pch = 16) 
lines(lowess(pisa$ESCS, pisa$MATH), col = 2, lwd = 2)
abline(lm(MATH ~ ESCS, pisa), col = 4, lwd = 2)


## ----scatter2-----------------------------------------------------------------------------------------------------------
plot(READ ~ ESCS, pisa, pch = 16)
lines(lowess(pisa$ESCS, pisa$READ), col = 2, lwd = 2)
abline(lm(READ ~ ESCS, pisa), col = 4, lwd = 2)
plot(MATH ~ READ, pisa, pch = 16)
lines(lowess(pisa$READ, pisa$MATH), col = 2, lwd = 2)
abline(lm(MATH ~ READ, pisa), col = 4, lwd = 2)


## ----xyplot1, fig.asp = 0.9---------------------------------------------------------------------------------------------
xyplot(MATH ~ ESCS | SCHOOLID,  data = pisa,  main = "MATH",
       panel = function(x, y){ 
         panel.xyplot(x, y)
         panel.loess(x, y)
         panel.lmline(x, y, lty = 2) } )


## ----xyplot2, fig.asp = 0.9---------------------------------------------------------------------------------------------
xyplot(READ  ~ MATH | SCHOOLID,  data = pisa, 
       main = "READ vs MATH",
       panel=function(x, y){ 
         panel.xyplot(x, y)
         panel.loess(x, y)
         panel.lmline(x, y, lty = 2) } )


## ----schooltypes, fig.asp = 0.7-----------------------------------------------------------------------------------------
par(mfrow = c(1, 2), pty = "s", mar = rep(3, 4))
boxplot(MATH ~ LYCEUM * TECH, pisa, 
        names = c("Other", "Lyceum", "Tech", ""), 
        xlim = c(0.5, 3.5), xlab = "School type")
boxplot(READ ~ LYCEUM * TECH, pisa, 
        names = c("Other", "Lyceum", "Tech", ""), 
        xlim = c(0.5, 3.5),  xlab = "School type")


## ----xyplot3, fig.asp = 0.9---------------------------------------------------------------------------------------------
xyplot(MATH ~ ESCS | SCHOOLID, data = subset(pisa, LYCEUM == TRUE), 
       main = "LYCEUM",
       panel = function(x, y){ 
         panel.xyplot(x, y)
         panel.loess(x, y)
         panel.lmline(x, y, lty = 2) } )
xyplot(MATH ~ ESCS | SCHOOLID, data = subset(pisa, LYCEUM == FALSE), 
       main = "NOT LYCEUM",
        panel = function(x, y){ 
          panel.xyplot(x, y)
          panel.loess(x, y)
          panel.lmline(x, y, lty = 2) } )


## ----margin1------------------------------------------------------------------------------------------------------------
MATH_schoolmean <- with(pisa, tapply(MATH, SCHOOLID, mean))
ESCS_schoolmean <- with(pisa, tapply(ESCS, SCHOOLID, mean))
plot(MATH_schoolmean ~ ESCS_schoolmean)
lines(lowess(ESCS_schoolmean, MATH_schoolmean), col = 2, lwd = 2)


## ----margin2------------------------------------------------------------------------------------------------------------
READ_schoolmean <- with(pisa, tapply(READ, SCHOOLID, mean))
plot(READ_schoolmean ~ ESCS_schoolmean)
lines(lowess(ESCS_schoolmean, READ_schoolmean), col = 2, lwd = 2)


## ----model1, fig.asp = 0.9----------------------------------------------------------------------------------------------
math.lm <- lm(MATH ~ ESCS, data = pisa)
summary(math.lm)
par(mfrow = c(2, 3), pty = "s", mar = rep(3, 4))
plot(math.lm, which = 1:6)


## ----groupres-----------------------------------------------------------------------------------------------------------
bwplot(getGroups(pisa) ~ resid(math.lm))
math.lm <- update(math.lm, . ~ . + factor(SCHOOLID))
bwplot(getGroups(pisa) ~ resid(math.lm))
anova(math.lm)


## ----anova--------------------------------------------------------------------------------------------------------------
math.lm <- update(math.lm, . ~ . + factor(SCHOOLID) * ESCS)
anova(math.lm) 


## ----lmlist-------------------------------------------------------------------------------------------------------------
pisa_nas <- colSums(is.na(pisa))
math.list <- lmList(MATH ~ ESCS | SCHOOLID, pisa[,  names(pisa)[pisa_nas==0]])
plot(intervals(math.list))


## ----multi1-------------------------------------------------------------------------------------------------------------
math.lme.null <- lme(MATH ~ 1, random = ~ 1 | SCHOOLID, data = pisa)
summary(math.lme.null)


## ----shrink-------------------------------------------------------------------------------------------------------------
math.list <- lmList(MATH ~ 1 | SCHOOLID, pisa[,  names(pisa)[pisa_nas==0]])
comp.math <- compareFits(coef(math.list), coef(math.lme.null))
plot(comp.math, mark = fixef(math.lme.null))


## ----multi2-------------------------------------------------------------------------------------------------------------
math.lme <- lme(MATH ~ ESCS + KIND + SEX + GRADE, random = ~ 1 | SCHOOLID, data = pisa) 
summary(math.lme)


## ----ML-----------------------------------------------------------------------------------------------------------------
math.lme.ML <-  lme(MATH ~ ESCS + KIND + SEX + GRADE, random = ~ 1 | SCHOOLID, 
                    data = pisa, method = "ML")
summary(math.lme.ML)


## ----CI-----------------------------------------------------------------------------------------------------------------
intervals(math.lme)


## ----components---------------------------------------------------------------------------------------------------------
fixef(math.lme)
ranef(math.lme)
coef(math.lme)
VarCorr(math.lme)


## ----res1---------------------------------------------------------------------------------------------------------------
r1 <- resid(math.lme, level = 1, type = "p")
bwplot(getGroups(math.lme) ~ r1)   


## ----res0---------------------------------------------------------------------------------------------------------------
r0 <- resid(math.lme, level = 0, type = "p")
bwplot(getGroups(math.lme) ~ r0)  


## ----diagn--------------------------------------------------------------------------------------------------------------
plot(math.lme, MATH ~ fitted(.))
plot(math.lme, resid(.) ~ ESCS, abline = 0)
plot(math.lme, resid(., type = "p") ~ fitted(.))
qqnorm(math.lme, ~ resid(.))  
qqnorm(math.lme, ~ resid(.) | SEX)  
qqnorm(math.lme, ~ ranef(.) | LYCEUM)


## ----tests--------------------------------------------------------------------------------------------------------------
V <- as.matrix(math.lme$apVar)
sigma <- exp(attr(math.lme$apVar, "Pars")[1])
se.lsigma <- sqrt(V[1, 1]) 
se.sigma <- se.lsigma * sigma
1 - pnorm(sigma / se.sigma)  


## ----LRTs---------------------------------------------------------------------------------------------------------------
math.slope <- update(math.lme, random = ~ ESCS | SCHOOLID)
w.oss <- 2 * (logLik(math.slope, REML = TRUE) - logLik(math.lme, REML = TRUE)) 
ogg.sim <- simulate.lme(math.lme, data = pisa, nsim = 10,
                      seed = 1988, m2 = list(random = ~ ESCS | SCHOOLID))
w.sim <- 2 * (ogg.sim$alt$REML[, 2] - ogg.sim$null$REML[, 2])


## ----marg---------------------------------------------------------------------------------------------------------------
math.marg <- gls(MATH ~ ESCS + KIND + SEX + GRADE ,data = pisa,
                 method = "REML", correlation = corCompSymm(form = ~ 1 | SCHOOLID))
cbind(math.marg$coef, fixef(math.lme))


## ----biv----------------------------------------------------------------------------------------------------------------
pisa2 <- rbind(pisa, pisa)
pisa2$Y <- c(pisa$MATH, pisa$READ)
pisa2$STUD <- c(1:nrow(pisa), 1:nrow(pisa))
pisa2$TYPOR <- factor(c(rep("math", nrow(pisa)), rep("read", nrow(pisa))))
mod.biv <-lme(Y ~ ESCS * TYPOR - 1 - ESCS, weights = varIdent(form = ~ 1 | TYPOR), 
             data = pisa2, correlation = corCompSymm(form = ~ 1 | SCHOOLID / STUD),
             random = ~ TYPOR - 1 | SCHOOLID)
summary(mod.biv)
intervals(mod.biv)


## ----betwith------------------------------------------------------------------------------------------------------------
math.lme <-  lme(MATH ~ ESCS + SEX + GRADE, random = ~ 1 | SCHOOLID, 
                 data = pisa, method = "ML")
math.lme2 <- lme(MATH  ~ CESCS + SEX + GRADE + ESCSM, random = ~ 1 | SCHOOLID, 
                 data = pisa, method = "ML")
mysummary <- function(object){
               assignInNamespace("print.correlation", 
                  function(x, title) return(), ns="nlme")  ###suppress correlation
              summary(object) 
              }
mysummary(math.lme2)


## ----anova2-------------------------------------------------------------------------------------------------------------
anova(math.lme, math.lme2) 


## ----update-------------------------------------------------------------------------------------------------------------
math.lme <- math.lme2


## ----studlevel----------------------------------------------------------------------------------------------------------
math.lme2 <- update(math.lme, .~. + KIND + FAMSTRUC + CHANGE + EARLY)
mysummary(math.lme2)


## ----schoolevel---------------------------------------------------------------------------------------------------------
math.lme2 <- update(math.lme, . ~ . + LYCEUM + TECH)
mysummary(math.lme2) 
math.lme <- math.lme2


## ----context------------------------------------------------------------------------------------------------------------
math.lme2 <- update(math.lme, .~. + SCHLSIZE + ISTMED + ISTLAR + PCGMED + PCGALT)
options(digits = 5)
mysummary(math.lme2)
math.lme2 <- update(math.lme, . ~ . + SCHLSIZE)
math.lme <- math.lme2


## ----schooldescr--------------------------------------------------------------------------------------------------------
math.lme2 <- update(math.lme, .~. + SCMATEDU + SCMATBUI + MEDDISC, na.action = na.omit)
options(digits = 5)
mysummary(math.lme2)
math.lme2 <- update(math.lme, .~. + MEDDISC, na.action = na.omit)


## ----final--------------------------------------------------------------------------------------------------------------
options(digits = 5)
mysummary(math.lme2)


## ----lme4---------------------------------------------------------------------------------------------------------------
library(lme4) 
data(Penicillin)
?Penicillin


## ----cross--------------------------------------------------------------------------------------------------------------
head(Penicillin)
xtabs( ~ sample + plate, Penicillin)


## ----lmer---------------------------------------------------------------------------------------------------------------
mod1 <- lmer(diameter ~ 1 + (1 | sample), data = Penicillin)
mod2 <- lmer(diameter ~ 1 + (1 | plate), data = Penicillin)
mod.vc <- lmer(diameter ~ 1 + (1 | sample) + (1 | plate), data = Penicillin)
summary(mod.vc)


## ----multi--------------------------------------------------------------------------------------------------------------
npos <- c(11, 16, 14, 2, 6, 1, 1, 4, 10, 22, 7, 1, 0, 0, 1, 6)
ntot <- c(36, 20, 19, 16, 17, 11, 5, 6, 37, 32, 19, 17, 12, 10, 9, 7)
treatment <- c(rep(1,8), rep(0,8))
clinic <-c(seq(8), seq(8))
booth <- data.frame(succ = npos, den = ntot, treat = treatment, cli = clinic)
booth$treat <- factor(booth$treat)
bh.glm0 <- glm(cbind(succ, den - succ) ~ treat, data = booth, family = binomial)
summary(bh.glm0)


## ----resplots, fig.asp = 0.95-------------------------------------------------------------------------------------------
library(boot)
glm.diag.plots(bh.glm0)


## ---- fig.asp = 0.95----------------------------------------------------------------------------------------------------
bh.glm1 <- update(bh.glm0, . ~ . + factor(cli))
glm.diag.plots(bh.glm1)


## ----AICbh--------------------------------------------------------------------------------------------------------------
AIC(bh.glm0, bh.glm1) 


## ----PQL----------------------------------------------------------------------------------------------------------------
bh.pql <- glmmPQL(cbind(succ, den - succ) ~ treat, 
                  random = ~ 1 | cli, 
                  data = booth, family = "binomial")
summary(bh.pql)
intervals(bh.pql)


## ----shrinklog----------------------------------------------------------------------------------------------------------
bh.glm1 <- glm(cbind(succ, den - succ) ~ factor(cli)- 1 + treat, data = booth,
              family = binomial)
tabe <- cbind(coef(bh.pql)[, 1], bh.glm1$coef[-9])
colnames(tabe) <- c("Random", "Fixed")
knitr::kable(tabe)


## ----AGH----------------------------------------------------------------------------------------------------------------
library(glmmML)
bh.ML1 <- glmmML(cbind(succ, den - succ) ~ treat, 
                 data = booth, 
                 family = binomial(logit),
                 cluster = booth$cli, method = "Laplace")
bh.ML10 <- glmmML(cbind(succ, den - succ) ~ treat, 
                  data = booth, 
                  family = binomial(logit),
                  cluster = booth$cli, method = "ghq", n.points = 10) 
bh.ML1
bh.ML10


## ----Cauchy-------------------------------------------------------------------------------------------------------------
bh.cauchy.ML <- glmmML(cbind(succ, den - succ) ~ treat, 
                       data = booth, 
                       family = binomial(logit),
                       cluster = booth$cli, method = "Laplace", prior = "cauchy")
bh.cauchy.ML 


## ----glmer--------------------------------------------------------------------------------------------------------------
library(lme4)
bh.ML4.10 <- glmer(cbind(succ, den - succ) ~ treat + (1 | cli), 
                   data = booth,
                   family = binomial(logit), nAGQ = 10)
bh.ML4.10


## ----probit-------------------------------------------------------------------------------------------------------------
bh.ML4.10 <- glmer(cbind(succ, den - succ) ~ treat + (1 | cli), data = booth,
                   family = binomial(probit), nAGQ = 10)
bh.ML4.10


## ----glmmTMB, warning = FALSE-------------------------------------------------------------------------------------------
library(glmmTMB)
booth$cli <- factor(booth$cli)
bh.TMB <- glmmTMB(succ / den ~ I(treat==1) + (1 | cli), 
                  data = booth, family = "binomial",
                  weights = den)
summary(bh.TMB)


## ----epildata-----------------------------------------------------------------------------------------------------------
data(epil2, package = "glmmTMB")
help(epil2)
epil2$subject <- factor(epil2$subject)


## ----visuale, fig.asp = 0.9---------------------------------------------------------------------------------------------
xyplot(y ~ period |subject,  
       data = epil2, 
       panel = function(x, y){ 
         panel.xyplot(x, y) } )


## ----twomet, warning = FALSE, message = FALSE---------------------------------------------------------------------------
mod.pois.pql <- glmmPQL(y ~ Base * trt + Age + Visit, 
                        random = ~ Visit | subject,
                        data = epil2, family = "poisson")
mod.pois.glmer <- glmer(y ~ Base * trt + Age + Visit + (Visit | subject),
                        data = epil2, family = "poisson") 
mod.pois.pql
mod.pois.glmer


## ----negbin, warning=FALSE----------------------------------------------------------------------------------------------
mod.pois <- glmmTMB(y ~ Base * trt + Age + Visit + (Visit | subject), 
                    data = epil2, family = "poisson")   
summary(mod.pois)    
mod.nbin1 <- glmmTMB(y ~ Base * trt + Age + Visit + (Visit | subject),
                     data = epil2, family = "nbinom1")  
anova(mod.pois, mod.nbin1) 
summary(mod.nbin1) 


## ----guImmun------------------------------------------------------------------------------------------------------------
library(mlmRev)
data(guImmun)
head(guImmun)
summary(guImmun)


## ----guIglm-------------------------------------------------------------------------------------------------------------
guI.glm <- glm(immun ~ kid2p + mom25p + ord  +
                ethn + momEd + husEd +
                momWork + rural + pcInd81, 
                family = binomial,  data = guImmun) 


## ----guire1-------------------------------------------------------------------------------------------------------------
guI.pql.comm <- glmmPQL(immun ~ kid2p + mom25p + ord  + 
                        ethn + momEd + husEd +
                        momWork + rural + pcInd81,
                        family = binomial, data = guImmun, random = ~ 1 | comm)


## ----guire2-------------------------------------------------------------------------------------------------------------
guI.pql.mom <- glmmPQL(immun ~ kid2p + mom25p + ord  +
                       ethn + momEd + husEd +
                       momWork + rural + pcInd81,
                       family = binomial, data = guImmun, random = ~ 1 | mom)


## ----guire3-------------------------------------------------------------------------------------------------------------
guI.pql3 <- glmmPQL(immun ~ kid2p + mom25p + ord  + ethn + momEd + husEd +
                    momWork + rural + pcInd81,
                    family = binomial, data = guImmun, random = ~ 1 | comm / mom)


## ----guire4, warning = FALSE--------------------------------------------------------------------------------------------
guI.TMB <- glmmTMB(immun ~ kid2p + mom25p + ord  + ethn + momEd + husEd +
                   momWork + rural + pcInd81 + (1 | comm / mom),
                   data = guImmun, family = binomial(logit))


## ----compare------------------------------------------------------------------------------------------------------------
mat.coef <- cbind(coef(guI.glm), guI.pql.comm$coef$fixed, guI.pql.mom$coef$fixed, 
                 guI.pql3$coef$fixed, fixef(guI.TMB)$cond)
se.pql3 <- diag(guI.pql3$varFix)^.5   
se.pql.comm <- diag(guI.pql.comm$varFix)^.5  
se.pql.mom <- diag(guI.pql.mom$varFix)^.5  
se.glm <- diag(vcov(guI.glm))^.5
se.lap <- sqrt(diag(vcov(guI.TMB)$cond))
mat.se <- cbind(se.glm, se.pql.comm, se.pql.mom, se.pql3, se.lap)
mat <- mat.coef / mat.se
options(digits = 2)
colnames(mat) <- c("glm", "PQL-comm", "PQL-family", "PQL-3 levels", "MLE-3 levels")
knitr::kable(mat)


## ----VC-----------------------------------------------------------------------------------------------------------------
VarCorr(guI.pql.comm)
VarCorr(guI.pql.mom)
VarCorr(guI.pql3)
VarCorr(guI.TMB)
qqnorm(guI.pql3, ~ ranef(., level = 2)) #family
qqnorm(guI.pql3, ~ ranef(., level = 1)) #community  


## ----comp, warning = FALSE----------------------------------------------------------------------------------------------
guI.ML.mom <- glmmML(immun ~ kid2p + mom25p + ord  +
                     ethn + momEd + husEd +
                     momWork + rural + pcInd81,
                     cluster = guImmun$mom,
                     method = "ghq", n.points = 10, 
                     data = guImmun, family = binomial(logit))
guI.ML.mom2 <- glmer(immun ~ kid2p + mom25p + ord  + ethn + momEd + husEd +
                     momWork + rural + pcInd81 + (1 | mom), nAGQ = 10,
                     data = guImmun, family = binomial(logit))


## ----final table--------------------------------------------------------------------------------------------------------
mat2 <- cbind(mat[, 3], guI.ML.mom$coeffi /  guI.ML.mom$coef.sd, 
             fixef(guI.ML.mom2) / sqrt(diag(vcov(guI.ML.mom2)))) 
colnames(mat2)<- c("PQL-family", "MLE-family", "MLE2-family")
knitr::kable(mat2)

