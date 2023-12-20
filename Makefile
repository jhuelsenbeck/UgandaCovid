LIBS     =   -lm

CFLAGS   =  -O4 -std=gnu++20 -Ofast

CC       =   g++

OBJECTS  =   main.o CondLikeJob.o CondLikeJobMngr.o Container.o MathCache.o Mcmc.o McmcInfo.o MetaData.o Model.o Msg.o Node.o Probability.o RandomVariable.o RateMatrix.o Threads.o TransitionProbabilities.o TransitionProbabilitiesMngr.o Tree.o UserSettings.o

PROGS    = cov

all:		$(PROGS)

cov:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o cov
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Container.o:	Container.cpp
		$(CC) $(CFLAGS) -c Container.cpp

CondLikeJob.o:	CondLikeJob.cpp
		$(CC) $(CFLAGS) -c CondLikeJob.cpp

CondLikeJobMngr.o:	CondLikeJobMngr.cpp
		$(CC) $(CFLAGS) -c CondLikeJobMngr.cpp

MathCache.o:	MathCache.cpp
		$(CC) $(CFLAGS) -c MathCache.cpp

Mcmc.o:	Mcmc.cpp
		$(CC) $(CFLAGS) -c Mcmc.cpp

McmcInfo.o:	McmcInfo.cpp
		$(CC) $(CFLAGS) -c McmcInfo.cpp

MetaData.o:	MetaData.cpp
		$(CC) $(CFLAGS) -c MetaData.cpp

Model.o:	Model.cpp
		$(CC) $(CFLAGS) -c Model.cpp

Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Node.o:	Node.cpp
		$(CC) $(CFLAGS) -c Node.cpp

Probability.o:	Probability.cpp
		$(CC) $(CFLAGS) -c Probability.cpp

RandomVariable.o:	RandomVariable.cpp
		$(CC) $(CFLAGS) -c RandomVariable.cpp

RateMatrix.o:	RateMatrix.cpp
		$(CC) $(CFLAGS) -c RateMatrix.cpp

Threads.o:	Threads.cpp
		$(CC) $(CFLAGS) -c Threads.cpp

TransitionProbabilities.o:	TransitionProbabilities.cpp
		$(CC) $(CFLAGS) -c TransitionProbabilities.cpp

TransitionProbabilitiesMngr.o:	TransitionProbabilitiesMngr.cpp
		$(CC) $(CFLAGS) -c TransitionProbabilitiesMngr.cpp

Tree.o:	Tree.cpp
		$(CC) $(CFLAGS) -c Tree.cpp

UserSettings.o:	UserSettings.cpp
		$(CC) $(CFLAGS) -c UserSettings.cpp

clean:		
		rm -f *.o
