Latest_version = "0.0.1.4"

# Install AVG R Package
Avg_RPackage_path<-file.path(args[5],"AvG_R_Package", paste0("AvantGardeDIA_", Latest_version,".tar.gz"))

## Install AvantGarde Package if it is not already installed
instaled_pckgs<-installed.packages()
instaled_pckgs<-instaled_pckgs[,1]
if (!is.element("AvantGardeDIA",instaled_pckgs) || packageVersion("AvantGardeDIA") < Latest_version ){
  install.packages(Avg_RPackage_path,repos = NULL)
}
