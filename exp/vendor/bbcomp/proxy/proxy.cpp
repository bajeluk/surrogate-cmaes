
// black box library header and function (pointer) declarations
#include <bbcomplib.h>
DLL_DECLARATIONS

// additional headers
#include "socket.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#ifdef _WIN32
	#pragma comment(lib, "ws2_32.lib")
#else
	#include <sys/select.h>
#endif


#define DEFAULT_PORT 7000


using namespace std;


string i2s(long value)
{
	stringstream ss;
	ss << value;
	return ss.str();
}

string d2s(double value, int precision = 20)
{
	stringstream ss;
	ss.precision(precision);
	ss << value;
	return ss.str();
}

string array2s(vector<double> const& value, int precision = 20)
{
	stringstream ss;
	ss.precision(precision);
	ss << "[" << value[0];
	for (size_t i=1; i<value.size(); i++) ss << " " << value[i];
	ss << "]";
	return ss.str();
}

int parseInt(const char* input, size_t& length)
{
	char* endp;
	int ret = strtol(input, &endp, 10);
	length = endp - input;
	return ret;
}

double parseDouble(const char* input, size_t& length)
{
	char* endp;
	double ret = strtod(input, &endp);
	length = endp - input;
	return ret;
}

vector<double> parseVector(const char* input, size_t& length)
{
	length = 0;
	vector<double> ret;
	while (input[length] == ' ') length++;
	if (input[length] != '[') { length = 0; return ret; }
	length++;
	while (input[length] == ' ') length++;
	if (input[length] == ']') { length++; return ret; }
	while (true)
	{
		char* endp;
		double value = strtod(input + length, &endp);
		if (endp == input + length) { length = 0; ret.clear(); return ret; }
		length = endp - input;
		ret.push_back(value);
		while (input[length] == ' ') length++;
		if (input[length] == ']') { length++; return ret; }
	}
}

string parseString(const char* input, size_t& length)
{
	length = 0;
	while (input[length] == ' ') length++;
	if (input[length] != '\"') { length = 0; return ""; }
	length++;
	size_t start = length;
	while (input[length] != '\"' && input[length] != 0) length++;
	if (input[length] != '\"') { length = 0; return ""; }
	length++;
	return string(input + start, length - start - 1);
}

bool parseComma(const char* input, size_t& length)
{
	length = 0;
	while (input[length] == ' ') length++;
	if (input[length] == ',') { length++; return true; }
	else return false;
}

bool parseClosingBracket(const char* input, size_t& length)
{
	length = 0;
	while (input[length] == ' ') length++;
	if (input[length] != ')') return false;
	length++;
	while (input[length] == ' ') length++;
	return (input[length] == 0);
}

// parse line into the form
//   functionname ( parameter , parameter , ... )
// where parameter is a number, an array of doubles
// in square brackets, or a double quoted string
string processLine(string line)
{
	// find opening bracket and function name
	size_t ob = line.find('(');
	if (ob == string::npos) return "ERROR,\"opening bracket not found\"";
	string functionname = line.substr(0, ob);
	while (! functionname.empty() && functionname[functionname.size() - 1] == ' ') functionname.erase(functionname.size() - 1);
	size_t pos = ob + 1;

	// case distinction
	size_t length;
	if (functionname == "configure")
	{
		// parse parameters
		int history = parseInt(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'history'\"";
		pos += length;
		if (! parseComma(line.c_str() + pos, length)) return "ERROR,\"comma expected between parameters\"";
		pos += length;
		string logfilepath = parseString(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'logfilepath'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = configure(history, logfilepath.c_str());

		// return result
		return i2s(result);
	}
	else if (functionname == "login")
	{
		// parse parameters
		string username = parseString(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'username'\"";
		pos += length;
		if (! parseComma(line.c_str() + pos, length)) return "ERROR,\"comma expected between parameters\"";
		pos += length;
		string password = parseString(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'password'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = login(username.c_str(), password.c_str());

		// return result
		return i2s(result);
	}
	else if (functionname == "numberOfTracks")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = numberOfTracks();

		// return result
		return i2s(result);
	}
	else if (functionname == "trackName")
	{
		// parse parameters
		int trackindex = parseInt(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'trackindex'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		const char* result = trackName(trackindex);

		// return result
		return string(result);
	}
	else if (functionname == "setTrack")
	{
		// parse parameters
		string trackname = parseString(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'trackname'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = setTrack(trackname.c_str());

		// return result
		return i2s(result);
	}
	else if (functionname == "numberOfProblems")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = numberOfProblems();

		// return result
		return i2s(result);
	}
	else if (functionname == "setProblem")
	{
		// parse parameters
		int problemID = parseInt(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'problemID'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = setProblem(problemID);

		// return result
		return i2s(result);
	}
	else if (functionname == "dimension")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = dimension();

		// return result
		return i2s(result);
	}
	else if (functionname == "numberOfObjectives")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = numberOfObjectives();

		// return result
		return i2s(result);
	}
	else if (functionname == "budget")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = budget();

		// return result
		return i2s(result);
	}
	else if (functionname == "evaluations")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int result = evaluations();

		// return result
		return i2s(result);
	}
	else if (functionname == "evaluate")
	{
		// parse parameters
		vector<double> point = parseVector(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'point'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int nObj = numberOfObjectives();
		if (nObj == 0) return "0";
		vector<double> value(nObj);
		int result = evaluate(&point[0], &value[0]);

		// return result
		if (result == 0) return "0";
		if (nObj == 1) return i2s(result) + "," + d2s(value[0]);
		else return i2s(result) + "," + array2s(value);
	}
	else if (functionname == "history")
	{
		// parse parameters
		int index = parseInt(line.c_str() + pos, length);
		if (length == 0) return "ERROR,\"failed to parse parameter 'index'\"";
		pos += length;
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		int dim = dimension();
		if (dim == 0) return "0";
		int nObj = numberOfObjectives();
		if (nObj == 0) return "0";
		vector<double> point(dim);
		vector<double> value(nObj);
		int result = history(index, &point[0], &value[0]);

		// return result
		if (result == 0) return "0";
		string ret = i2s(result);
		ret += ",";
		ret += array2s(point);
		ret += ",";
		if (nObj == 1) ret += d2s(value[0]);
		else ret += array2s(value);
		return ret;
	}
	else if (functionname == "errorMessage")
	{
		// parse parameters
		if (! parseClosingBracket(line.c_str() + pos, length)) return "ERROR,\"closing bracket expected at end of command\"";

		// call library
		const char* result = errorMessage();

		// return result
		return string(result);
	}
	else
	{
		return "ERROR,\"unknown function name\"";
	}
}

int main(int argc, char** argv)
{
	// connect to black box library
	DLL_LOAD
	if (DLL_LOAD_STATUS < 0)
	{
		cout << "error loading library" << endl;
		exit(EXIT_FAILURE);
	}

	// process command line arguments
	unsigned short port = DEFAULT_PORT;
	if (argc > 1)
	{
		long p = strtol(argv[1], NULL, 10);
		if (p >= 0 && p <= 65535) port = p; 
	}
	cout << "proxy ready on port " << port << endl;

	// wait for connections and requests
	Socket listener(port, "127.0.0.1");
	Socket* sock = NULL;
	string readbuffer;
	while (true)
	{
		fd_set socks;
		FD_ZERO(&socks);
		listener.populateSet(socks);
		if (sock != NULL) sock->populateSet(socks);
		timeval timeout;
		timeout.tv_sec = 1;
		timeout.tv_usec = 0;

		// select call
		int count = select(FD_SETSIZE, &socks, NULL, NULL, &timeout);

		if (count < 0)
		{
			listener.close();
			break;
		}
		else if (count > 0)
		{
			if (listener.isInSet(socks))
			{
				if (sock == NULL)
				{
					// accept the connection
					sock = listener.accept();
				}
				else
				{
					// refuse further connections (since the black box library is stateful)
					Socket* socket = listener.accept();
					socket->close();
					delete socket;
				}
			}
			else if (sock != NULL && sock->isInSet(socks))
			{
				// read package
				char buffer[4096];
				size_t n = sock->read(buffer, 4096);
				if (n == 0)
				{
					delete sock;
					sock = NULL;
				}
				else
				{
					size_t start = readbuffer.size();
					readbuffer.append(buffer, n);
					size_t pos = readbuffer.find('\n', start);
					if (pos != string::npos)
					{
						// process line (one function call)
						string line = readbuffer.substr(0, pos);
						readbuffer.erase(0, pos + 1);
						string reply = processLine(line);
						if (! reply.empty())
						{
							// send reply
							reply += "\n";
							size_t sent = 0;
							while (sent < reply.size())
							{
								size_t n = sock->write(reply.c_str() + sent, reply.size() - sent);
								if (n == 0)
								{
									sock->close();
									delete sock;
									break;
								}
								sent += n;
							}
						}
					}
				}
			}
		}
	}

	// quit cleanly
	DLL_CLOSE
}
