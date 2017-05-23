
#pragma once

#include <string>

#ifdef _WIN32
	#define _WINSOCK_DEPRECATED_NO_WARNINGS
	#include <ws2tcpip.h> 
	#include <windows.h>
#else
	#include <sys/types.h>
#endif


//
// Generic encapsulation of a stream (TCP/IP) socket.
// The socket supports
//  1. listening (server),
//  2. (synchronous or asynchronous) communication (client and server).
// It should work on all major operating systems (Linux/Mac/Windows).
//
class Socket
{
public:
	// uninitialized (unconnected) socket
	Socket();

	// connect to a remote host
	Socket(std::string remoteUrl, unsigned short remotePort);

	// listen on a local port, possibly restricted to certain IP (e.g. "127.0.0.1")
	explicit Socket(unsigned short localPort, std::string remoteIP = "");

	// close the connection
	~Socket();


	// connect to remote host
	bool connect(std::string remoteUrl, unsigned short remotePort);

	// start listening
	bool listen(unsigned short localPort, std::string remoteIP = "");

	// close the connection (makes isGood() return false)
	void close();


	// check whether the connection is fine
	inline bool isGood() const
	{ return (m_handle > 0); }

	// to whom is this socket connected?
	bool getPeer(unsigned int& ip, unsigned short& port);

	// check whether the socket has data ready for an asynchronous read operation
	bool hasData(unsigned int timeoutSeconds = 0, unsigned int timeoutUSeconds = 0);

	// asynchronous read, may return size zero if no data is pending
	std::size_t read(char* buffer, std::size_t buffersize);

	// elementary write operations are always synchronous
	// (this is how the OS implements it)
	std::size_t write(const char* buffer, std::size_t buffersize);

	// Accept a pending connection (if hasData() is false then this call is blocking)
	// and return a newly allocated socket object representing this connection.
	// The caller is responsible for deleting the returned socket object.
	Socket* accept();

	// add the socket's descriptor to the set
	void populateSet(fd_set& set) const;

	// check whether the socket is listed in the set
	bool isInSet(fd_set& set) const;

private:
	// used only to disambiguate constructors
	enum disambiguate
	{
		constructFromHandle,
	};

	// construction from OS socket handle
	explicit Socket(disambiguate dis, int handle);

	// POSIX socket handle
	int m_handle;

	// network library initialization
	static void initialize();
	static bool m_initialized;
};
