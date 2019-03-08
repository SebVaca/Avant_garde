Latest_version = "0.0.1.0"

# Install AVG R Package
Avg_RPackage_path<-file.path(args[5],"AvG_R_Package", paste0("AvantGardeDIA_", Latest_version,".tar.gz"))

## Install AvantGarde Package if it is not already installed
a<-installed.packages()
packages<-a[,1]
if (!is.element("AvantGardeDIA",packages) || packageVersion("AvantGardeDIA") < Latest_version ){
  install.packages(Avg_RPackage_path,repos = NULL)
}
