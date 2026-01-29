require(roxygen2)
require(devtools)
require(pkgdown)

# Document
setwd("~/GitHub/NICHES")
document()
check()

# Check installation from local
setwd("~/GitHub")
install('NICHES')

# Update website
setwd("~/Documents/GitHub/NICHES")
pkgdown::build_site()
pkgdown::build_article(name = "07 Spatiotemporal NICHES",lazy = T)

# Install from release ?
install_github('msraredon/NICHES', 
               ref = 'v1.2.4',
               force = T)
