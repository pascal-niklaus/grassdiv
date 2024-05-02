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

options(digits = 3)
close_all()
rm(list = ls())

d0 <- rbind(
    data.frame(hv = 1, read.csv("derived_data/CE_shoots_hv_1.csv")),
    data.frame(hv = 2, read.csv("derived_data/CE_shoots_hv_2.csv")),
    data.frame(hv = 3, read.csv("derived_data/CE_shoots_hv_3.csv")))

## the pyrenee soil has a productivity that is approx 3 times smaller.
## adjust for this using an empirical factor

d0$k <- ifelse(d0$soil == "P", 3, 1)
d0 <- transform(
    d0,
    acccomp = sorted.code(d0$tacc,d0$cacc,collapse=":"),
    spcomp = sorted.code(d0$tsp,d0$csp,collapse=":"),
    cshoot = cshoot * k,
    tshoot = tshoot * k,
    shoots = shoots * k
)

fromUX <- function(x) {
    x0 <- par("usr")[1]
    x1 <- par("usr")[2]
    dx <- x1 - x0
    (x - x0) / dx
}




######################################################################
###
### monocultures only

d <- subset(d0,  tacc == cacc)

for (h in 1:3) {
    simple.heading(sprintf("hv = %d", h), pre=1, post=1)

    cat("---------- using pot as replicate\n")
    d.aov <- aov(log(shoots) ~ tsp * tsite * soil,
                 subset = (h == hv),
                 data = d)
    print(summary(d.aov))


    cat("---------- using species as replicate\n")
    d.aov <- aov(log(shoots) ~ tsp + tsite * soil + Error(tsp:tsite:soil),
                 subset = (h == hv),
                 data = d)
    print(summary(d.aov))
}

## sum of all three harvests
d123 <- aggr(d, c("pot","tacc","soil","tsite","tsp"),"shoots=mean(shoots)")
d123.aov <- aov(log(shoots) ~ tsp * tsite * soil,
                data = d123)
print(summary(d123.aov))

######################################################################
###
### Figure 1a - shoot biomass

d0 <- rbind(
    data.frame(hv = 1, read.csv("derived_data/CE_shoots_hv_1.csv")),
    data.frame(hv = 2, read.csv("derived_data/CE_shoots_hv_2.csv")),
    data.frame(hv = 3, read.csv("derived_data/CE_shoots_hv_3.csv")))

d2 <- transform(
    d0,
    shoot = cshoot + tshoot,
    comp = sorted.code(d0$tsp, d0$csp, collapse = "·", unique = TRUE)
)
d2 <- aggr(d2, c("pot", "soil", "spdiv", "comp"), c("shoot=sum(shoot)"))

d2.aggr <- aggr(d2, c("soil", "spdiv", "comp"), c("shoot=mean(shoot)", "se.shoot=se(shoot)"))
d2.aggr <- transform(
    d2.aggr,
    colr = ifelse(spdiv == 1, "white", " gray")
)

cairo_pdf("Fig1a.pdf", width = 6, height = 6)
#
par(mar = c(8, 4, .5, .2))
xpos <- barplot(
    d2.aggr$shoot,
    space = groupSpace(d2.aggr, c("soil", "spdiv")),
    ylim = c(0, 14),
    col = d2.aggr$colr,
    las = 2,
    ylab = "",
    cex.lab =1.4
)
xy.errbar(xpos, d2.aggr$shoot, yplus = d2.aggr$se.shoot, add = TRUE)
mtext(expression("Shoot biomass (g pot"^{-1}*")"), side = 2, cex = 1.2, line = 2)
text(xpos, rep(-.5, length(xpos)), d2.aggr$comp, srt=90, adj = c(1, 0.5), offset = 0.5, xpd = NA)
segments(x0 = xpos[1]-.3, x1 = xpos[10]+.3, y0=-3.5, xpd=NA)
segments(x0 = xpos[11] - .3, x1 = xpos[20] + .3, y0 = -3.5, xpd = NA)
text(x = mean(range(xpos[1:10])), y = -4.2, "Alps", xpd=NA)
text(x = mean(range(xpos[11:20])), y = -4.2, "Pyrenees", xpd=NA)
dev.off()

######################################################################
###
### Figure 1b - net diversity effects

d2 <- transform(
    d0,
    shoot = cshoot + tshoot,
    comp = sorted.code(d0$tsp, d0$csp, collapse = "·", unique = TRUE),
    acomp = sorted.code(d0$tacc, d0$cacc, collapse = "·", unique = TRUE)
)

d2 <- aggr(d2, c("pot", "soil", "spdiv", "tacc", "cacc", "comp", "acomp"), c("shoot=sum(shoot)"))

d2.monos <- subset(d2, tacc == cacc)
d2.monos <- aggr(d2.monos, c("soil", "acomp=tacc"), "shoot=mean(shoot)")

d2.mix <- subset(d2, spdiv == 2)

key <- paste(d2.mix$acomp,d2.mix$soil,sep="·")
d2.mix$monomean <- unlist(lapply(
    strsplit(key, "·", fixed = TRUE),
    function(x) {
        tmp <- subset(d2.monos, soil == x[3])
        x <- x[1:2]
        mean(tmp$shoot[tmp$acomp %in% x])
    }
))
d2.mix <- transform(d2.mix, NE = shoot - monomean)
d2.NE <- aggr(d2.mix, c("soil", "comp"), c("NE=mean(NE)", "se.NE=se(NE)"))

cairo_pdf("Fig1b.pdf", width = 6, height = 6)
##
par(mar = c(8, 4, .5, .2))
xpos <- barplot(
    d2.NE$NE,
    space = groupSpace(d2.NE, "soil"),
    ylim = c(-0.5, 2.5),
    col = "gray",
    las = 2,
    ylab = "",
    cex.lab =1.4
)
xy.errbar(
    xpos, d2.NE$NE,
    yplus = (d2.NE$NE > 0) * (d2.NE$se.NE),
    yminus = (d2.NE$NE < 0) * (d2.NE$se.NE),
    add = TRUE
)
mtext(expression("Net diversity effect (g pot"^{-1}*")"), side = 2, cex = 1.2, line = 2)
text(xpos, rep(-.5, length(xpos)), d2.NE$comp, srt=90, adj = c(1, 0.5), offset = 0.5, xpd = NA)
segments(x0 = xpos[1]-.3, x1 = xpos[6]+.3, y0=-1, xpd=NA)
segments(x0 = xpos[7] - .3, x1 = xpos[12] + .3, y0 = -1, xpd = NA)
text(x = mean(range(xpos[1:6])), y = -1.2, "Alps", xpd=NA)
text(x = mean(range(xpos[7:12])), y = -1.2, "Pyrenees", xpd=NA)
##
dev.off()

######################################################################
###
### new combined figure

library(gridBase)
library(grid)

cairo_pdf("Fig1.pdf", width = 10, height = 6)

## empty top figure without margin
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NULL, xlim = 0:1, ylim = 0:1)
## push inner figure viewport
vps <- gridBase::baseViewports()
pushViewport(vps$inner)

## split available area into grid
lay <- grid.layout(
    nrow = 4, ncol = 5,
    widths = unit(
        c(4, 2.2, 5.5, 1.3, 0.5),
        c("lines", "null", "lines", "null", "lines")
    ),
    heights <- unit(
        c(1, 1, 5, 1),
        c("lines", "null", "lines", "lines")
    )
)

pushViewport(viewport(layout = lay))

vpl <- vpList(
    viewport(layout.pos.col = 2, layout.pos.row = 2, name = "div"),
    viewport(layout.pos.col = 4, layout.pos.row = 2, name = "NE"),
    viewport(layout.pos.col = 2, layout.pos.row = 3, name = "leg.div"),
    viewport(layout.pos.col = 4, layout.pos.row = 3, name = "leg.NE")
)
pushViewport(vpl)

## ----------------------- original data plot ------------------------

seekViewport("div")

par(plt = gridPLT())
xpos <- barplot(
    d2.aggr$shoot,
    space = groupSpace(d2.aggr, c("soil", "spdiv")),
    ylim = c(0, 14),
    col = d2.aggr$colr,
    las = 2,
    ylab = "",
    xaxs="r",
    cex.lab = 1.4
)
xy.errbar(xpos, d2.aggr$shoot, yplus = d2.aggr$se.shoot, add = TRUE)
mtext(expression("Shoot biomass (g pot"^{-1}*")"), side = 2, cex = 1.2, line = 2)

grid.text(d2.aggr$comp,x=fromUX(xpos), y=unit(-.5, "lines"), just="right", rot=70)

grid.lines(x = fromUX(range(xpos[1:10]) + c(-.5, .5)), y = unit(-4.8, "char"), gp = gpar(lwd = 2))
grid.lines(x = fromUX(range(xpos[11:20]) + c(-.5, .5)), y = unit(-4.8, "char"), gp = gpar(lwd = 2))

grid.text("Alps", x = fromUX(mean(range(xpos[1:10]))), y = unit(-5.8, "char"))
grid.text("Pyrenees", x = fromUX(mean(range(xpos[11:20]))), y = unit(-5.8, "char"))


## ------------------------- net effect plot -------------------------

repl1st <- function(x, repl) {
    c(repl, x[-1])
}


seekViewport("NE")
par(plt = gridPLT(), new=TRUE)

xpos <- barplot(
    d2.NE$NE,
    space = repl1st(groupSpace(d2.NE, "soil"),1),
    ylim = c(-0.5, 2.5),
    col = "gray",
    las = 2,
    ylab = "",
    xaxs="r",
    cex.lab = 1.4
)

xy.errbar(
    xpos, d2.NE$NE,
    yplus = (d2.NE$NE > 0) * (d2.NE$se.NE),
    yminus = (d2.NE$NE < 0) * (d2.NE$se.NE),
    cap = 0.03,
    add = TRUE
)
mtext(expression("Net diversity effect (g pot"^{-1}*")"), side = 2, cex = 1.2, line = 2)

grid.text(d2.NE$comp,x=fromUX(xpos), y=unit(-.5, "lines"), just="right", rot=70)

grid.lines(x = fromUX(range(xpos[1:6]) + c(-.5, .5)), y = unit(-4.8, "char"), gp = gpar(lwd = 2))
grid.lines(x = fromUX(range(xpos[7:12]) + c(-.5, .5)), y = unit(-4.8, "char"), gp = gpar(lwd = 2))

grid.text("Alps", x = fromUX(mean(range(xpos[1:7]))), y = unit(-5.8, "char"))
grid.text("Pyrenees", x = fromUX(mean(range(xpos[7:12]))), y = unit(-5.8, "char"))


## --------------------------------- ---------------------------------
popViewport(0)

dev.off()


######################################################################
###
### Additive partitioning

d1 <- aggr(
    d0,
    c(
        "pot", "soil",
        "tacc", "cacc",
        "tsite", "csite",
        "tsp", "csp",
        "spdiv", "accdiv",
        "spcomp", "acccomp"
    ),
    c("tshoot=sum(tshoot)", "cshoot=sum(cshoot)")
)

dorig <- transform(d1, sym = ifelse(tsite == csite, "sym", "allo"))
dorig <- restorefactors(dorig)

dorig <- stck(dorig,
    factors = c("pot", "soil", "spdiv", "spcomp", "accdiv", "acccomp", "sym"),
    to.stack = c(
        "shoot=cshoot,tshoot",
        "acc=tacc,cacc",
        "sp=tsp,csp"
    ),
    cat.names = c("plant=target,competitor")
)

dorig<- dorig[order(dorig$pot, dorig$plant), ]
rownames(dorig) <- NULL

dorig$acccomp <- sorted.code(dorig$acccomp,split = ":", collapse = ":", unique = TRUE)

addpart(shoot ~ acccomp / acc + pot, groups = ~sym, data = dorig)
