corr <- function(directory, threshold = 0) {
        source("complete.R")
        com <- complete(directory)
        
        # subset the data.frame according to the threshold for the nobs
        data <- com[com$nobs > threshold, ]
        
        # result is a numeric vector
        result <- numeric(0)
        
        # for each data point, read CSV, calculate the cor and append to the result
        for(i in data$id) {
                csv<-read.csv(paste(getwd(), "/", directory, "/", sprintf("%03d", i), ".csv", sep=''))
                
                # logical vector of valid rows
                tf <- !is.na(csv$sulfate) & !is.na(csv$nitrate)
                
                # subset of valid rows
                x <- csv[tf, ]  
                
                result <- c(result, cor(x$sulfate, x$nitrate))
        }
        
        result
}