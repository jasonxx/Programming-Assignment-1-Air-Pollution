complete<-function(directory, id=1:332){
       nobs = numeric(0)
        for(i in id) {
                csv<-read.csv(paste(getwd(), "/", directory, "/", sprintf("%03d", i), ".csv", sep=''))
                nobs = c(nobs, sum(!is.na(csv$sulfate) & !is.na(csv$nitrate)))
        }
        data.frame(id=id, nobs=nobs)
}