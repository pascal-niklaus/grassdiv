######################################################################
###
### Local adaptation
###
### Pascal A. Niklaus
###
### History:
### - 16-02-2022 file created

library(pascal)
library(emmeans)
library(gridBase)
library(grid)

options(digits=3)
close_all()
rm(list=ls())

usr2dev <- function(x) {
    w <- par("pin")[1] / diff(par("usr")[1:2])
    h <- par("pin")[2] / diff(par("usr")[3:4])
    asp <- w / h
    angle <- 180 / pi * atan(x / asp)
    angle
}

sitelabel <- function(x) {
    ifelse(as.character(x) == "F", "Alps", "Pyrenees")
}

spname <- function(x) {
    unname(
        sapply(as.character(x), function(x) {
            if (x == "Des") {
                "Deschampsia"
            } else if (x == "Nar") {
                "Nardus"
            } else if (x == "Phl") {
                "Phleum"
            } else if (x == "Poa") {
                "Poa"
            } else {
                stop("Unknown species name", x, "!")
            }
        })
    )
}

d0 <- rbind(
    data.frame(hv=1, read.csv("derived_data/CE_shoots_hv_1.csv")),
    data.frame(hv=2, read.csv("derived_data/CE_shoots_hv_2.csv")),
    data.frame(hv=3, read.csv("derived_data/CE_shoots_hv_3.csv")))
d0 <- restorefactors(d0)

## the pyrenee soil has a productivity that is approx 3 times smaller.
## adjust for this using an empirical factor

d0$k <- ifelse(d0$soil == "P", 3, 1)

d0 <- transform(d0,
                cshoot = cshoot * k,
                tshoot = tshoot * k,
                shoots = shoots * k)

d <- aggr(d0,
          c("pot", "tacc", "cacc", "soil", "tsite", "csite", "tsp",
            "csp", "spdiv", "accdiv"),
          c("shoots=sum(shoots)", "tshoot=sum(tshoot)", "cshoot=sum(cshoot)"))

## Give monos and mixtures the same weight

da <- aggr(
    d,
    c("tsp", "csp", "tsite", "soil"),
    c("y=mean(log(tshoot))", "n=length(tshoot)")
)

## LA analysis:

## - terms for target and competitor accession origin and soil origin.
## - For LA test to co-occurring plant species, monos excluded
## - log-transformed cum. shoot biomass of target,
## - target sp, comp sp, target acc origin, comp acc orig,
##   target acc x origin accession origin as fixed terms;
## - interaction accessions and soil as error term
## - i.e. 48 species—site origin-combinations as replicates fortests.

## - For LA to soil env, sp mono kept.
## - trg_sp, cmp_sp, torig, soil, torig x sorig
## - 64 sp—site origin combinations as replicates

## - for LA analysis by species: no tsp
## - random term with 4 competitor species x 4 site—origin combs


######################################################################
###
### LA to soil

## three possibilities:
## - use tacc x  soil as replicate -> hopeless, because too few reps
## - use tacc x competition species as replicate
if (FALSE) {
    m1 <- aov(
        log(tshoot) ~ tsp + csp + tsite * soil + Error(tacc:soil),
        data = d
    )
    summary(m1)
    ## Error: tacc:soil
    ##            Df Sum Sq Mean Sq F value  Pr(>F)
    ## tsp         3  17.09    5.70   47.58 7.6e-06 ***
    ## tsite       1   0.37    0.37    3.10    0.11
    ## soil        1   0.04    0.04    0.33    0.58
    ## tsite:soil  1   0.25    0.25    2.07    0.18
    ## Residuals   9   1.08    0.12

    ## Error: Within
    ##            Df Sum Sq Mean Sq F value Pr(>F)
    ## csp         3   3.61   1.203    16.2  3e-09 ***
    ## Residuals 157  11.63   0.074

    m2 <- aov(
        log(tshoot) ~ tsp + csp + tsite * soil + Error(tacc:csp:soil),
        data = d
    )
    summary(m2)
    ## Error: tacc:csp:soil
    ##            Df Sum Sq Mean Sq F value  Pr(>F)
    ## tsp         3  17.09    5.70   88.43 < 2e-16 ***
    ## csp         3   3.61    1.20   18.68 1.9e-08 ***
    ## tsite       1   0.37    0.37    5.76   0.020 *
    ## soil        1   0.04    0.04    0.62   0.435
    ## tsite:soil  1   0.25    0.25    3.85   0.055 .
    ## Residuals  54   3.48    0.06

    ## Error: Within
    ##            Df Sum Sq Mean Sq F value Pr(>F)
    ## Residuals 112   9.23  0.0824


    m3 <- aov(
        log(tshoot) ~ tsp + cacc + tsite * soil + Error(tacc:cacc:soil),
        data = d
    )
    summary(m3)
    ## Error: tacc:cacc:soil
    ##             Df Sum Sq Mean Sq F value  Pr(>F)
    ## tsp          3  17.09    5.70   84.83 < 2e-16 ***
    ## cacc         7   5.55    0.79   11.80 3.1e-11 ***
    ## tsite        1   0.59    0.59    8.81  0.0037 **
    ## soil         1   0.04    0.04    0.59  0.4427
    ## tsite:soil   1   0.25    0.25    3.69  0.0571 .
    ## Residuals  114   7.65    0.07

    ## Error: Within
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Residuals 48   2.89  0.0603
}

m2a <- aov(
    y ~ tsp + csp + tsite * soil,
    data = da
)
print(summary(m2a))


    ##             Df Sum Sq Mean Sq F value  Pr(>F)
    ## tsp          3   7.65   2.551  108.45 < 2e-16 ***
    ## csp          3   1.54   0.513   21.82 2.2e-09 ***
    ## tsite        1   0.24   0.240   10.19  0.0024 **
    ## soil         1   0.02   0.019    0.80  0.3755
    ## tsite:soil   1   0.17   0.166    7.05  0.0104 *
    ## Residuals   54   1.27   0.024


cairo_pdf("Fig2.pdf", width = 10, height = 6)

xex <- .3

par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NULL, xlim = 0:1, ylim = 0:1)
vps <- gridBase::baseViewports()
pushViewport(vps$inner)

## split available area into grid
lay <- grid.layout(
    nrow = 4, ncol = 8,
    widths = unit(
        c(5, 1, 2, 1, 1, 1, 1, 0.5),
        c("lines", "null", "lines", "null", "null", "null", "null", "lines")
    ),
    heights <- unit(
        c(3, 1, 5, 1),
        c("lines", "null", "lines", "lines")
    )
)

pushViewport(viewport(layout = lay))

vpl <- vpList(
    viewport(layout.pos.col = 2, layout.pos.row = 2, name = "all"),
    viewport(layout.pos.col = 4, layout.pos.row = 2, name = "Des"),
    viewport(layout.pos.col = 5, layout.pos.row = 2, name = "Nar"),
    viewport(layout.pos.col = 6, layout.pos.row = 2, name = "Phl"),
    viewport(layout.pos.col = 7, layout.pos.row = 2, name = "Poa"),
    viewport(layout.pos.col = 2:7, layout.pos.row = 3, name = "xlegend"),
    viewport(layout.pos.col = 1, layout.pos.row = 2, name = "ylegend")
)
pushViewport(vpl)

seekViewport("all")

par(plt = gridPLT())

daa <- aggr(da, c("soil", "tsite"), c("y=mean(y)", "se=se(y)"))
yvec <- seq(2, 9, by = 1)
daa <- as.data.frame(emmeans(m2a, specs = c("tsite", "soil")))
daa <- transform(
    daa,
    x = as.numeric(soil) + scale(as.numeric(tsite), scale = FALSE) * .05
)

plot(
    emmean ~ x,
    data = daa,
    pch = ifelse(daa$tsite == "F", 16, 1),
    cex = 2,
    yaxt = "n", xaxt = "n", yaxs = "i",
    xlab = "", ylab = "",
    xlim = range(daa$x) + c(-xex,+xex),
    ylim = log(range(yvec))
)
# mtext(expression("Shoot biomass (g)"), side = 2, cex = 1.2, line = 3)
# mtext(expression("Soil origin"), side = 1, cex = 1.2, line = 3)
mtext("All species", side = 3, cex = 1.2, line = 1)
xy.errbar(daa$x, daa$emmean, yerr = daa$SE, add = TRUE)
for (org in daa$tsite) {
    tmp <- subset(daa, tsite == org)
    lines(tmp$x, tmp$emmean)
}
left.axis(at = log(yvec), labels = sprintf("%.1f", yvec), las = 2)
right.axis(at = log(yvec), labels = rep("", length(yvec)), las = 2, tck = .025)
bottom.axis(at = as.numeric(daa$soil), labels = sitelabel(daa$soil), las = 1)
tmp <- aggr(daa, "tsite", c("x=mean(x)", "y=mean(emmean)"))
for (si in suc(daa$tsite)) {
    tmp <- subset(daa, tsite == si)
    slope <- unname(coef(lm(tmp$emmean ~ tmp$x))[2])
    text(mean(tmp$x),
        mean(tmp$emmean),
        sitelabel(tmp$tsite)[1],
        cex = 1.0, pos = 3,
        srt = usr2dev(slope)
    )
}
pvalixn <- summary(m2a)[[1]]["tsite:soil", "Pr(>F)"]
mtext(sigStars(pvalixn), side = 3, line = -1.5)

text(1.5, log(6.5), "Accession origin",cex=1)


## plot by species
for (sp in suc(da$tsp)) {
    seekViewport(sp)
    par(plt = gridPLT(), new = TRUE)

    tmp <- subset(da, tsp == sp)
    tmp.aov <- aov(y ~ csp + tsite * soil, data = tmp)
    pvalixn <- summary(tmp.aov)[[1]]["tsite:soil","Pr(>F)"]

    tmpa <- as.data.frame(emmeans(tmp.aov, specs = c("tsite", "soil")))
    tmpa <- transform(
        tmpa,
        x = as.numeric(soil) + scale(as.numeric(tsite), scale = FALSE) * .05,
        y = emmean
    )

    plot(
        y ~ x,
        data = tmpa,
        pch = ifelse(tmpa$tsite == "F", 16, 1),
        cex = 2,
        yaxt = "n", xaxt = "n", yaxs = "i",
        xlab = "", ylab = "",
        xlim = range(tmpa$x) + c(-xex,+xex),
        ylim = log(range(yvec))
    )
    mtext(spname(sp), side = 3, cex = 1.2, line = 1)
    bottom.axis(at = as.numeric(tmpa$soil), labels = sitelabel(tmpa$soil), las = 1)
    left.axis(at = log(yvec), labels = rep("", length(yvec)), las = 2, tck = .025)
    right.axis(at = log(yvec), labels = rep("", length(yvec)), las = 2, tck = .025)

    mtext(sigStars(pvalixn), side = 3, line = -1.5)

    print(tmpa)
    xy.errbar(tmpa$x, tmpa$emmean, yerr = tmpa$SE, add = TRUE)
    for (org in tmpa$tsite) {
        tmp <- subset(tmpa, tsite == org)
        lines(tmp$x, tmp$emmean)
    }
}

seekViewport("xlegend")
grid.text("Soil origin", x = 0.5, y = 0.3, gp=gpar(fontsize=18))
seekViewport("ylegend")
grid.text("Shoot biomass (g)", x = 0.3, y = 0.5, rot=90,gp=gpar(fontsize=18))

popViewport(0)
dev.off()


######################################################################
###
### LA to
