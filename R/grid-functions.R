arrangeSideBySide <- function(object1, object2){
	grid.newpage()
	lvp <- viewport(x=0,
			y=0.05,
			width=unit(0.50, "npc"),
			height=unit(0.95, "npc"), just=c("left", "bottom"),
			name="lvp")
	pushViewport(lvp)
	nfigs1 <- length(object1$condlevels[[1]])
	nfigs2 <- length(object2$condlevels[[1]])
	stopifnot(length(nfigs1) == length(nfigs2))
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	object1$layout <- c(1, nfigs1)
	print(object1, newpage=FALSE, prefix="plot1", more=TRUE)
	upViewport(0)
	lvp2 <- viewport(x=0.5,
			 y=0.05,
			 width=unit(0.50, "npc"),
			 height=unit(0.95, "npc"), just=c("left", "bottom"),
			 name="lvp2")
	pushViewport(lvp2)
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	object2$layout <- c(1, nfigs1)
	object2
	print(object2, newpage=FALSE, prefix="plot2", more=TRUE)
}
