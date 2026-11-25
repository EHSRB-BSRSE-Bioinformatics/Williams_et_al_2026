source("R/revision/gsea_unique_ids.R")
gsea_environment()
capture.output(dput(formals(fgsea::fgseaMultilevel)),file=file.path(gsea_docs,"fgsea_release_defaults.txt"))
capture.output(sessionInfo(),file=file.path(gsea_docs,"fgsea_environment_sessionInfo.txt"))
namespaces <- sort(loadedNamespaces())
root <- normalizePath(getwd())
source("R/provenance.R")
packages <- do.call(rbind,lapply(namespaces,function(p)data.frame(package=p,
  version=as.character(getNamespaceVersion(p)),
  path={
    location <- if(p=="base") file.path(R.home("library"),"base") else getNamespaceInfo(asNamespace(p),"path")
    relative <- archive_relative_path(location, root)
    if (identical(relative, location)) basename(location) else relative
  })))
write.table(packages,file.path(gsea_docs,"loaded_package_versions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(.libPaths(),file.path(gsea_docs,"library_paths.txt"))
cat("Verified fgsea",as.character(getNamespaceVersion("fgsea")),"on",R.version.string,"\n")
