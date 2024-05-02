######################################################################
###
### Biodiversity-Effects
###
### Pascal A. Niklaus
###
### History:
### - 16-02-2022 file created

library(pascal)
library(pdiv)
library(emmeans)
library(gridBase)
library(grid)

options(digits = 3)
close_all()
rm(list = ls())

dev_from_usr_x <- function(x) {
    x0 <- par("usr")[1]
    x1 <- par("usr")[2]
    dx <- x1 - x0
    (x - x0) / dx
}

angle_from_usr_slope <- function(x) {
    w <- par("pin")[1] / diff(par("usr")[1:2])
    h <- par("pin")[2] / diff(par("usr")[3:4])
    asp <- w / h
    angle <- 180 / pi * atan(x / asp)
    angle
}

d0 <- rbind(
    data.frame(hv = 1, read.csv("derived_data/CE_shoots_hv_1.csv")),
    data.frame(hv = 2, read.csv("derived_data/CE_shoots_hv_2.csv")),
    data.frame(hv = 3, read.csv("derived_data/CE_shoots_hv_3.csv")))

## the pyrenee soil has a productivity that is approx 3 times smaller.
## adjust for this using an empirical factor

d0$k <- ifelse(d0$soil == "P", 3, 1)

d0 <- transform(d0,
                cshoot = cshoot * k,
                tshoot = tshoot * k,
                shoots = shoots * k)

d0 <- aggr(
    d0, c("pot", "tacc", "cacc", "soil"),
    c("cshoot=sum(cshoot)", "tshoot=sum(tshoot)")
)


d <- transform(d0,
    shoot = cshoot + tshoot,
    acccomp = sorted.code(tacc, cacc, collapse = "|", unique = TRUE)
)

d2 <- stck(
    d, c("pot", "soil", "acccomp"),
    to.stack = c("acc=cacc,tacc", "shoot=cshoot,tshoot")
)

## take sum of C and T in monoculture
d2 <- aggr(d2, c("pot", "soil", "acccomp", "acc"), "shoot=sum(shoot)")

d2.ap <- addpart(shoot ~ acccomp / acc + pot, groups = ~soil, data = d2)

d2.ap <- transform(
    d2.ap,
    acc1 = sapply(strsplit(acccomp, "|", fixed = TRUE), function(x) x[1]),
    acc2 = sapply(strsplit(acccomp, "|", fixed = TRUE), function(x) c(x, x)[2])
)

d2.ap <- transform(
    d2.ap,
    site1 = substr(acc1, 5, 6),
    site2 = substr(acc2, 5, 6),
    sp1 = substr(acc1, 1, 3),
    sp2 = substr(acc2, 1, 3)
)


# exclude accession monocultures
d2.ap <- subset(d2.ap, acc1 != acc2)

d2.ap <- transform(
    d2.ap,
    accsym = ifelse(site1 == site2, "sym", "allo"),
    soilsym = ifelse(site1 == site2 & site1 == soil, "2sym",
        ifelse(site1 == soil | site2 == soil, "1sym", "0sym")
    ),
    spmono = ifelse(sp1 == sp2, "mono", "mix"),
    acc1 = NULL,
    acc2 = NULL
)
usr_angle_to_dev <- function(x) {
    w <- par("pin")[1] / diff(par("usr")[1:2])
    h <- par("pin")[2] / diff(par("usr")[3:4])
    asp <- w / h
    angle <- 180 / pi * atan(x / asp)
    angle
}



d2.ap <- moveColumns(d2.ap, c("sp1", "sp2", "spmono", "site1", "site2","accsym","soilsym"), after="soil")

d2.ap <- aggr(
    d2.ap,
    c("soil", "sp1", "sp2", "spmono",
      "site1", "site2", "accsym", "soilsym", "acccomp"),
    c("CE=mean(CE.shoot)", "SE=mean(SE.shoot)")
)

######################################################################
###
### Plot CE, SE
###
### 1) by allo/sym in species monos
### 2) by sym relative to soil
### 3) by allo/sym in species mixtures

######################################################################
###
### Soil environment: CE/SE vs # of plants on home soil

cairo_pdf("Fig3.pdf", width = 10, height = 6)

par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NULL, xlim = 0:1, ylim = 0:1)
vps <- gridBase::baseViewports()
pushViewport(vps$inner)

## split available area into grid
lay <- grid.layout(
    nrow = 4, ncol = 7,
    widths = unit(
        c(6, 1, 1, 1, 4, .7, 3),
        c("lines", "null", "lines", "null", "lines", "null", "lines")
    ),
    heights <- unit(
        c(3, 1, 5, 1),
        c("lines", "null", "lines", "lines")
    )
)

pushViewport(viewport(layout = lay))

vpl <- vpList(
    viewport(layout.pos.col = 2, layout.pos.row = 2, name = "soil"),
    viewport(layout.pos.col = 4, layout.pos.row = 2, name = "mix"),
    viewport(layout.pos.col = 6, layout.pos.row = 2, name = "mono"),
    viewport(layout.pos.col = 2, layout.pos.row = 3, name = "xlegendsoil"),
    viewport(layout.pos.col = 4, layout.pos.row = 3, name = "xlegendmix"),
    viewport(layout.pos.col = 6, layout.pos.row = 3, name = "xlegendmono"),
    viewport(layout.pos.col = 2:4, layout.pos.row = 1, name = "titlemix"),
    viewport(layout.pos.col = 6, layout.pos.row = 1, name = "titlemono"),
    viewport(layout.pos.col = 1, layout.pos.row = 2, name = "ylegend")
)
pushViewport(vpl)

####################################################################
###
### species mixtures vs number of accessions on home soil (0, 1, 2)

seekViewport("soil")
par(plt = gridPLT())

d <- subset(d2.ap, spmono == "mix")

d_CE_soil <- aggr(
    d,
    c("soilsym"),
    c("CE=mean(CE)", "se.CE=se(CE)", "SE=mean(SE)", "se.SE=se(SE)")
)

d_CE_soil <- transform(d_CE_soil,
    nhome = as.numeric(substring(soilsym, 1, 1))
)

shft <- -.03
xext <- .2
with(
    d_CE_soil,
    {
        plot(
            y = CE, x = nhome + shft,
            xlim = range(nhome) + c(-xext, +xext),
            ylim = c(0, 1),
            xaxs = "i", yaxs = "i",
            xaxt = "n", yaxt = "n",
            xlab = "", ylab = "",
            pch = 16
        )
        xy.errbar(
            y = CE, x = nhome + shft, yplus = se.CE,
            pch = 16, add = TRUE, cap = .07
        )
        lines(y = CE, x = nhome + shft, pch = 1, type = "b")
        points(y = SE, x = nhome - shft, pch = 1)
        lines(y = SE, x = nhome - shft, pch = 1, type = "b")
        xy.errbar(
            y = SE, x = nhome - shft, yminus = se.SE,
            pch = 16, add = TRUE, cap = .07
        )
    }
)
left.axis()
right.axis(labels = FALSE, tck = .025)

tmp <- d_CE_soil[-1, ]
rotCE <- angle_from_usr_slope(
    coef(lm(CE ~ nhome, data = tmp))["nhome"]
)
rotSE <- angle_from_usr_slope(
    coef(lm(SE ~ nhome, data = tmp))["nhome"]
)

text(1.5, mean(tmp$CE), "CE", pos = 3, srt = rotCE)
text(1.5, mean(tmp$SE), "SE", pos = 1, srt = rotSE)

bottom.axis(at = 0:2)

seekViewport("xlegendsoil")
grid.text(
    "Number of accessions\non soil of origin",
    x = 0.5, y = 0.3, gp = gpar(fontsize = 15)
)


seekViewport("titlemix")
grid.text(
    "Species mixtures",
    x = 0.5, y = 0.5, gp = gpar(fontsize = 18)
)

## mtext("Species mixtures", side=3, line = .5,cex = 1.2)

####################################################################
###
### species mixtures vs sym/allo accessions

seekViewport("mix")
par(plt = gridPLT(), new = TRUE)

d <- subset(d2.ap, spmono != "mono")
d_CE_mix <- aggr(
    d,
    c("spmono", "accsym"),
    c("CE=mean(CE)", "se.CE=se(CE)", "SE=mean(SE)", "se.SE=se(SE)")
)
d_CE_mix <- transform(d_CE_mix, orig_div = ifelse(accsym == "sym", 0, 1))
with(
    d_CE_mix,
    {
        plot(y = CE, x = orig_div + shft,
             xlim = range(orig_div) + c(-xext, +xext),
             ylim = c(0, 1),
             pch = 16,
             xaxs = "i", yaxs = "i",
             xaxt = "n", yaxt = "n",
             xlab = "", ylab = "")
        abline(h = 0)
        xy.errbar(
            y = CE, x = orig_div + shft,
            yplus = se.CE, cap = .07,
            pch = 16, add = TRUE
        )
        lines(y = CE, x =  orig_div + shft, pch = 1, type = "b")
        points(y = SE, x = orig_div - shft, pch = 1)
        lines(y = SE, x = orig_div - shft, pch = 1, type = "b")
        xy.errbar(
            y = SE, x = orig_div - shft,
            yminus = se.SE,
            pch = 16, cap = .07, add = TRUE
        )
    }
)
left.axis(labels = FALSE, tck=.025)
right.axis(labels = FALSE, tck=.025)
bottom.axis(at = 0:1)

tmp <- d_CE_mix
rotCE <- angle_from_usr_slope(
    coef(lm(CE ~ orig_div, data = tmp))["orig_div"]
)
rotSE <- angle_from_usr_slope(
    coef(lm(SE ~ orig_div, data = tmp))["orig_div"]
)

text(.5, mean(tmp$CE), "CE", pos = 3, srt = rotCE)
text(.5, mean(tmp$SE), "SE", pos = 1, srt = rotSE)

seekViewport("xlegendmix")
grid.text(
    "Number of\naccessions origins",
    x = 0.5, y = 0.3, gp = gpar(fontsize = 15)
)


####################################################################
###
### accession monos vs mixtures in sp monos

seekViewport("mono")
par(plt = gridPLT(), new = TRUE)


shft <- -.02
d <- subset(d2.ap, spmono == "mono")
d_CE_mono <- aggr(
    d,
    c("spmono"),
    c("CE=mean(CE)", "se.CE=se(CE)", "SE=mean(SE)", "se.SE=se(SE)")
)

with(
    d_CE_mono,
    {
        plot(
            y = CE, x = shft,
            xlim = c(-xext, +xext), ylim = c(-.6, .6),
            xaxs = "i", yaxs = "i",
            xaxt = "n", yaxt = "n",
            xlab = "", ylab = "",
            pch = 16
        )
        abline(h = 0)
        xy.errbar(
            y = CE, x = shft, yerr = se.CE,
            pch = 16, add = TRUE, cap = 0.07
        )
        lines(y = CE, x =  shft, pch = 1, type = "b")
        points(y = SE, x = -shft, pch = 1)
        lines(y = SE, x = -shft, pch = 1, type = "b")
        xy.errbar(
            y = SE, x = -shft, yerr = se.SE,
            pch = 16, add = TRUE, cap = 0.07
        )
    }
)
left.axis()

text(shft, d_CE_mono$CE, "CE", adj = c(1.5,0.5))
text(-shft, d_CE_mono$SE, "SE", adj = c(-0.5,0.5))

seekViewport("xlegendmono")
grid.text(
    "Accession mixtures\n(relative to acc. monocultures)",
    x = 0.5, y = 0.3, gp = gpar(fontsize = 15)
)

seekViewport("titlemono")
grid.text(
    "Species monocultures",
    x = 0.5, y = 0.5, gp = gpar(fontsize = 18)
)


######################################################################
###
### Y legend

seekViewport("ylegend")
grid.text(
    expression("Component of net biodiversity effect"),
    x = unit(0.15,"npc"), y = 0.5, rot=90,gp=gpar(fontsize=18)
)
grid.text(
    expression("(g "*pot^{-1}*")"),
    x = unit(0.15,"npc") + unit(1, "lines"), y = 0.5, rot=90,gp=gpar(fontsize=18)
)

popViewport(0)
dev.off()
