# dataset.BindingDB

https://www.bindingdb.org/bind/chemsearch/marvin/SDFdownload.jsp?all_download=yes

### Import database into Oracle on Docker 

#### Pull docker 
```
docker pull sath89/oracle-12c
```

#### Run with 8080 and 1521 ports opened:
```
docker run --name bind -d -p 8080:8080 -p 1521:1521 sath89/oracle-12c

docker run --name bind -d -p 8080:8080 -p 1521:1521 -v `pwd`/oracle:/u01/app/oracle -v `pwd`/oracle_data:/oracle_data sath89/oracle-12c

docker exec -it bind bash 
```

#### Reset and Restart 
```
docker kill $(docker ps -a -q); docker rm $(docker ps -a -q); 
rm -rf oracle/*
docker volumne prune 
```  

#### Connect db with sqlplus  
```
cd oracle_data
gzip -d BindingDB_All_Oracle.dmp
sqlplus system/oracle
```

```
sqlplus system/oracle @reset_database.sql
impdp -parfile imp_all.txt ignore=y
```
