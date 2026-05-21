#include "global.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#ifdef _WIN32
#include <windows.h>
#else
#define MAX_PATH 4096
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#endif

double RandomDouble(double min, double max)
{
    return min + double(rand())/(double(RAND_MAX/(max - min)));
}

void RenameConsole(const string &title)
{
//	LPCWSTR t = title.c_str();
//    SetConsoleTitleA(title.c_str());
}

string CreateUniqueFileName(const string &filename, const string &ext)
{
    string name = filename + ext;

    for (int i = 1; ifstream(name).is_open(); ++i)
    {
        name = filename + '(' + to_string(i) + ')' + ext;
    }

    return name;
}

string CreateFolder(string &name)
{
    string num = "";
#ifdef _WIN32
    char curDir[MAX_PATH] = "";

    if (!GetCurrentDirectoryA(MAX_PATH, curDir))
    {
        cerr << "Error getting current directory: #" << GetLastError();
    }

    string basePath = string(curDir) + "\\" + name;
    string dirPath = basePath;

    for (int i = 1; !CreateDirectoryA(dirPath.c_str(), NULL); ++i)
    {
        num = "(" + to_string(i) + ")";
        dirPath = basePath + num;
    }
#else
    string dirName = name;

    for (int i = 1; mkdir(dirName.c_str(), 0755) != 0; ++i)
    {
        num = "(" + to_string(i) + ")";
        dirName = name + num;
    }
#endif
    return name + num;
}

string CreateDir(const string &name)
{
    string dirName = name;
#ifdef _WIN32
    string dirPath = name;

    for (int i = 1; !CreateDirectoryA(dirPath.c_str(), NULL); ++i)
    {
        dirPath = name + "(" + to_string(i) + ")";
    }

    dirName = dirPath + "/";
#else
    for (int i = 1; mkdir(dirName.c_str(), 0755) != 0; ++i)
    {
        dirName = name + "(" + to_string(i) + ")";
    }

    dirName += "/";
#endif
    return dirName;
}

void EraseConsoleLine(int lenght)
{
    cout << '\r';

    for (int i = 0; i < lenght; ++i)
    {
        cout << ' ';
    }

    cout << '\r';
}

double DegToRad(double deg)
{
    return (deg*M_PI)/180;
}

double RadToDeg(double rad)
{
    return (rad*180)/M_PI;
}
