.onLoad <- function (libname, pkgname)
{
    
    utils::globalVariables ("model.matrix")
    utils::globalVariables ("t.test")
    utils::globalVariables ("p.adjust")
    utils::globalVariables ("aov")
    utils::globalVariables ("as.formula")
    utils::globalVariables ("factors.x")
    utils::globalVariables ("factors.y")
    invisible ()

}