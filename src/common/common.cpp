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
    char curDir[MAX_PATH] = ""; // current directory
#ifdef _WIN32
    char newDir[MAX_PATH]; // created directory

    // get the current directory, and store it
    if (!GetCurrentDirectoryA(MAX_PATH, curDir))
    {
        cerr << "Error getting current directory: #" << GetLastError();
    }

    strcat(curDir, "\\");
    strcpy(newDir, curDir);
    strcat(newDir, name.c_str());

    char dirN[MAX_PATH];
    strcpy(dirN, newDir);
    string num = "";

    for (int i = 1; !CreateDirectoryA(dirN, NULL); ++i)
    {
        num = "(" + to_string(i) + ")";
        string dirName = newDir + num;
        strcpy(dirN, dirName.c_str());
    }

    name += num;
    return curDir;
#else
    // Linux: create directory, append (N) if already exists
    char cwd[MAX_PATH];
    if (!getcwd(cwd, MAX_PATH))
    {
        cerr << "Error getting current directory" << endl;
        return "";
    }

    string base = string(cwd) + "/" + name;
    string dirPath = base;
    string num = "";

    int rc = mkdir(dirPath.c_str(), 0755);
    if (rc != 0 && errno == EEXIST)
    {
        for (int i = 1; ; ++i)
        {
            num = "(" + to_string(i) + ")";
            dirPath = base + num;
            rc = mkdir(dirPath.c_str(), 0755);
            if (rc == 0 || errno != EEXIST) break;
        }
    }

    name += num;
    return string(cwd) + "/";
#endif
}

string CreateDir(const string &name)
{
    string dirName;
#ifdef _WIN32
    char cdir[MAX_PATH];
    cdir[0] = '\0';
    strcat(cdir, name.c_str());

    char dirN[MAX_PATH];
    strcpy(dirN, cdir);

    for (int i = 1; !CreateDirectoryA(cdir, NULL); ++i)
    {
        string dirName = dirN;
        dirName += "(" + to_string(i) + ")";
        strcpy(cdir, dirName.c_str());
    }

    strcat(cdir, "\\");
    dirName = cdir;
#else
    dirName = name;
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
