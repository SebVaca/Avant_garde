# Avant-garde R package version
a<-installed.packages()
packages<-a[,1]
if(is.element("AvantGardeDIA",packages)) {
  	msg = paste("AvantGardeDIA R package version: ", packageVersion("AvantGardeDIA"))
  } else{
  	msg = "AvantGardeDIA R package is not installed."
  }
  
# R version
R.Version()$version.string
print(msg)
