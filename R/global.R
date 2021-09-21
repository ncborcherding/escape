.onLoad <- function (libname, pkgname)
{
    
    utils::globalVariables ("model.matrix")
    utils::globalVariables ("t.test")
    utils::globalVariables ("p.adjust")
    utils::globalVariables ("aov")
    utils::globalVariables ("as.formula")
    utils::globalVariables ("factors.x")
    utils::globalVariables ("factors.y")
    utils::globalVariables ("slot")
    utils::globalVariables ("stats")
    utils::globalVariables ("glm")
    utils::globalVariables ("wilcox.test")
    invisible ()

}