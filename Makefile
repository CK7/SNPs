bin/SNPs : bin/SNPs.o bin/ScafSNPs.o bin/common.o bin/ReadMapping.o bin/ReadMappingReader.o bin/ReadSequence.o bin/String.o
	g++ bin/SNPs.o bin/ScafSNPs.o bin/common.o bin/ReadMapping.o bin/ReadMappingReader.o bin/ReadSequence.o bin/String.o -o bin/SNPs

bin/SNPs.o : source/SNPs.cpp
	g++ -c source/SNPs.cpp -o bin/SNPs.o

bin/ScafSNPs.o : source/ScafSNPs.cpp
	g++ -c source/ScafSNPs.cpp -o bin/ScafSNPs.o

bin/common.o : source/common/common.cpp
	g++ -c source/common/common.cpp -o bin/common.o

bin/ReadMapping.o : source/common/ReadMapping.cpp
	g++ -c source/common/ReadMapping.cpp -o bin/ReadMapping.o

bin/ReadMappingReader.o : source/common/ReadMappingReader.cpp
	g++ -c source/common/ReadMappingReader.cpp -o bin/ReadMappingReader.o

bin/String.o : source/common/String.cpp
	g++ -c source/common/String.cpp -o bin/String.o

bin/ReadSequence.o : source/common/ReadSequence.cpp
	g++ -c source/common/ReadSequence.cpp -o bin/ReadSequence.o

clean:
	rm -f bin/*.o
