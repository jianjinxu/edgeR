#  ZZZ.R

.onAttach <- function(libname, pkgname)
#	Add User's Guide to Windows menu
#	Gordon Smyth
#	14 Sep 2009.
{
	if(.Platform$OS.type=="windows" && .Platform$GUI=="Rgui" ) {
		winMenuAddItem("Vignettes","edgeR","shell.exec(system.file(\"doc\",\"edgeR.pdf\",package=\"edgeR\"))")
	}
}
