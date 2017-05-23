
#include "socket.h"

#include <cstring>
#include <stdexcept>
#include <cerrno>
#include <csignal>

#ifdef _WIN32
	typedef int socklen_t;

	#ifndef MSG_NOSIGNAL
		#define MSG_NOSIGNAL 0
	#endif

	#include <Ws2tcpip.h>

	#pragma comment(lib, "ws2_32.lib")
#else
	#include <netinet/tcp.h>
	#include <netinet/in.h>
	#include <arpa/inet.h>
	#include <sys/socket.h>
	#include <sys/time.h>
	#include <netdb.h>
	#include <unistd.h>
	#ifndef MSG_NOSIGNAL
		#define MSG_NOSIGNAL SO_NOSIGPIPE
	#endif
#endif


#define DISABLE_NAGLE


using namespace std;


// static
bool Socket::m_initialized = false;


Socket::Socket()
: m_handle(0)
{
	initialize();
}

Socket::Socket(string remoteUrl, unsigned short remotePort)
: m_handle(0)
{
	initialize();

	if (! connect(remoteUrl, remotePort)) throw runtime_error("[Socket::Socket] failed to connect");
}

Socket::Socket(unsigned short localPort, string remoteIP)
: m_handle(0)
{
	initialize();

	if (! listen(localPort, remoteIP)) throw runtime_error("[Socket::Socket] failed to listen");
}

Socket::Socket(disambiguate dis, int handle)
: m_handle(handle)
{
	initialize();

	(void)dis;
}

Socket::~Socket()
{
	close();
}

// static
void Socket::initialize()
{
	if (! m_initialized)
	{
#ifdef _WIN32
		WSADATA info;
		int result = WSAStartup(2 /* API version 2.0 */, &info);
		if (result != 0) throw runtime_error("[Socket::initialize] windows socket API startup failed");
#else
		signal(SIGPIPE, SIG_IGN);
#endif
		m_initialized = true;
	}
}


bool Socket::connect(string remoteUrl, unsigned short remotePort)
{
	// refuse if socket is already connected
	if (m_handle > 0) return false;

	m_handle = ::socket(AF_INET, SOCK_STREAM, 0);
	if (m_handle <= 0) return false;

#ifdef DISABLE_NAGLE
	int flag = 1;
	setsockopt(m_handle, IPPROTO_TCP, TCP_NODELAY, (char*)&flag, sizeof(flag));
#endif

	hostent *host = gethostbyname(remoteUrl.c_str());
	if (host == NULL)
	{
		close();
		return false;
	}

	sockaddr_in addr;
	addr.sin_family = AF_INET;
	addr.sin_port = htons(remotePort);
	addr.sin_addr.s_addr = ((in_addr *)host->h_addr)->s_addr;
	memset(&(addr.sin_zero), 0, 8);

	int ret = ::connect(m_handle, (sockaddr*) & addr, sizeof(sockaddr));
	if (ret < 0) close();

	return isGood();
}

bool Socket::listen(unsigned short localPort, string remoteIP)
{
	// refuse if socket is already connected
	if (m_handle > 0) return false;

	m_handle = ::socket(AF_INET, SOCK_STREAM, 0);
	if (m_handle <= 0) return false;

#ifdef _WIN32
	char yes = 1;
#else
	int yes = 1;
#endif
	if (setsockopt(m_handle, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1)
	{
		close();
		return false;
	}

	sockaddr_in addr;
	addr.sin_family = AF_INET;
	addr.sin_port = htons(localPort);
	if (remoteIP.empty()) addr.sin_addr.s_addr = INADDR_ANY;
	else
	{
#ifdef _WIN32
		if (InetPton(AF_INET,remoteIP.c_str(), &addr.sin_addr) == 0)
#else
		if (inet_aton(remoteIP.c_str(), &addr.sin_addr) == 0)
#endif
		{
			close();
			return false;
		}
	}
	memset(&(addr.sin_zero), 0, 8);

	if (::bind(m_handle, (sockaddr*)&addr, sizeof(sockaddr)) == -1)
	{
		close();
		return false;
	}

	if (::listen(m_handle, 5) == -1)
	{
		close();
		return false;
	}

	return true;
}

void Socket::close()
{
#ifdef _WIN32
	if (m_handle > 0) ::closesocket(m_handle);
#else
	if (m_handle > 0) ::close(m_handle);
#endif
	m_handle = 0;
}

bool Socket::getPeer(unsigned int& ip, unsigned short& port)
{
	sockaddr_in s;
	socklen_t len = sizeof(sockaddr_in);
	int err = getpeername(m_handle, (sockaddr*)(&s), &len);
	if (err == 0)
	{
		ip = ntohl(s.sin_addr.s_addr);
		port = ntohs(s.sin_port);
		return true;
	}
	else
	{
		close();
		return false;
	}
}

bool Socket::hasData(unsigned int timeoutSeconds, unsigned int timeoutUSeconds)
{
	if (m_handle == 0) return false;

	timeval timeout;
	timeout.tv_sec = timeoutSeconds;
	timeout.tv_usec = timeoutUSeconds;

	fd_set read_socks;
	FD_ZERO(&read_socks);
	FD_SET(m_handle, &read_socks);

	int count = select(m_handle + 1, &read_socks, NULL, NULL, &timeout);
	// bool ret = (FD_ISSET(m_handle, &read_socks));
	bool ret = (count == 1);

	FD_ZERO(&read_socks);
	if (count == -1)
	{
		close();
		return false;
	}

	return ret;
}

size_t Socket::read(char* buffer, size_t buffersize)
{
	int ret = recv(m_handle, buffer, buffersize, 0);
	if (ret <= 0)
	{
		ret = 0;
		close();
	}
	return ret;
}

size_t Socket::write(const char* buffer, size_t buffersize)
{
	if (buffersize == 0) return 0;
	int ret = send(m_handle, buffer, buffersize, MSG_NOSIGNAL);
	if (ret <= 0)
	{
		ret = 0;
		close();
	}
	return ret;
}

Socket* Socket::accept()
{
	int handle = ::accept(m_handle, NULL, NULL);
	if (handle < 0) return NULL;
	return new Socket(constructFromHandle, handle);
}

void Socket::populateSet(fd_set& set) const
{
	FD_SET(m_handle, &set);
}

bool Socket::isInSet(fd_set& set) const
{
	return (FD_ISSET(m_handle, &set) != 0);
}
