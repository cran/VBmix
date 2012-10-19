// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h> /* close */
#include <netdb.h> /* gethostbyname */
#define INVALID_SOCKET -1
#define SOCKET_ERROR -1
#define closesocket(s) close(s)
typedef int SOCKET;
typedef struct sockaddr_in SOCKADDR_IN;
typedef struct sockaddr SOCKADDR;
typedef struct in_addr IN_ADDR;

#define ADDRESS "127.0.0.1"
#define CRLF		"\r\n"
#define PORT	 	1979
#define MAX_CLIENTS 	1

#define BUF_SIZE	8196

int main(char argc, char** argv) {
   SOCKET sock = socket(AF_INET, SOCK_STREAM, 0);
   SOCKADDR_IN sin = { 0 };
   struct hostent *hostinfo;

   if(sock == INVALID_SOCKET) {
      printf("socket()\n");
      exit(-1);
   }

	// get address
	char * used_address;
	if(argc < 2) {
		used_address = ADDRESS;
	} else {
		used_address = argv[1];
	}


   hostinfo = gethostbyname(used_address);
   if (hostinfo == NULL) {
      printf("unknown host %s\n", used_address);
      exit(-1);
   }

	printf("connected to %s\n", used_address);

   sin.sin_addr = *(IN_ADDR *) hostinfo->h_addr;
   sin.sin_port = htons(PORT);
   sin.sin_family = AF_INET;

   if(connect(sock,(SOCKADDR *) &sin, sizeof(SOCKADDR)) == SOCKET_ERROR) {
      printf("connect()\n");
      exit(-1);
   }

	char buffer[BUF_SIZE];
	char *ptr = buffer;
	int nprint;
	nprint = sprintf(ptr, "flood");
	ptr += nprint;
	*ptr++ = 0;
	
	while(1) {
		send(sock, buffer, strlen(buffer), 0);
		printf("has flooded\n");
		sleep(2);
	}

	return(0);

}
