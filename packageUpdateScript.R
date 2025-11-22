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
               ref = 'v1.2.0',
               auth_token = 'f7ea5d8790fe721ac0c9d5ef115d04068b19ed6d',
               force = T)
