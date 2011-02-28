// Copyright (C) 2011 Pierrick Bruneau, see README for full notice

#include <stdio.h>
#include <stdlib.h>
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


#define CRLF		"\r\n"
#define PORT	 	1979
#define MAX_CLIENTS 	1

#define BUF_SIZE	8196

int main() {
	SOCKET sock = socket(AF_INET, SOCK_STREAM, 0);
   	SOCKADDR_IN sin = { 0 };
	int tr=1;
	int n = 0;
   	char buffer[BUF_SIZE];
    SOCKADDR_IN csin = { 0 };
    socklen_t sinsize = sizeof csin;

   	sin.sin_addr.s_addr = htonl(INADDR_ANY);
   	sin.sin_port = htons(PORT);
   	sin.sin_family = AF_INET;



	setsockopt(sock,SOL_SOCKET,SO_REUSEADDR,&tr,sizeof(int));
	bind(sock,(SOCKADDR *) &sin, sizeof sin);
	listen(sock, MAX_CLIENTS);   

	int csock = accept(sock, (SOCKADDR *)&csin, &sinsize);
    
	while(1) {

		n = recv(csock, buffer, BUF_SIZE - 1, 0);
		if(n != 0) {
   			buffer[n] = 0;
			printf("%s\n", buffer);
		} else {
			printf("transmission ended\n");
			break;
		}
	}
	
	return(0);
}
