pollutantmean <- function(directory, pollutant, id = 1:332) {
        data1<-NULL
        for (i in id){
                csv<-read.csv(paste(getwd(), "/", directory, "/", sprintf("%03d", i), ".csv", sep=''))
                data1<-rbind(data1,csv)
        }
        mean(data1[[pollutant]], na.rm = TRUE)        
}