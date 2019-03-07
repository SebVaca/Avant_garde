a<-installed.packages()
packages<-a[,1]
if(is.element("AvantGardeDIA",packages)) {
  print(paste("AvantGardeDIA R package version: ", packageVersion("AvantGardeDIA")))
  }
