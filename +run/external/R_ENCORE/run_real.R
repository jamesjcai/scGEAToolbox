source('readzzz.R')



# cell x gene
data=readzzz('matrix.txt','genelist.txt')

library("ENCORE")

result1=encore_step1(data,thread=12)
#2 run consensus clustering

result2=encore_step2(data=data,result=result1,method="dbscan",thread=1)

# result1=encore_step1(data,dens=2000,thread=2,sample="TRUE",sample_size=100)

ggplot(data=result2$temp,aes(x=tsne_1,y=tsne_2,color=clusters))+geom_point()

markers=find_markers(data, temp=result2$temp, thread=10)


# big data
result1=encore_step1(data,dens=2000,thread=2,sample="TRUE",sample_size=100)
result2=encore_step2(result=result1, ac=c(23,42), k=5, sample="TRUE", thread=4)
result=encore_big(data=data, result=result2,result1=result1,thread=2)
ggplot(data=result,aes(x=tsne_1,y=tsne_2,color=clusters))+geom_point()


