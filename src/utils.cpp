#include "../include/utils.h"


void isOut(int &i, int &j, int num)
{
    if(j >= num)
    {
        i++;
        j -= j;
    }
}

double string2double(string a)
{
    stringstream b;
    double num;
    b << a;
    b >> num;
    return num;
}

vector31d GetData(string path)
{
    ifstream in;
    in.open(path);
    if(!in.is_open())
    {
        cerr << "open file incorrectly!" << endl;
        exit(1);
    }
    cout << "reading all the data in the file..." << endl;
    vector31d data;
    int j = 0, count = 0;
    while(!in.eof())
    {
        
        string tmp;
        Vector3d vec;
        getline(in, tmp, ',');
        if(count==0)
        {
            count++;
            continue;
        }
        bool islast = false;
        if(tmp.find('\n') != string::npos)
        {
            tmp = tmp.substr(0, tmp.find('n') - 1);
            islast = true;
        }
        double num = string2double(tmp);
        cout << num << "        test" << endl;
        vec(j, 0) = num;
        j++;
        count++;
        if(j % 3 == 0)
        {
            data.push_back(vec);
            j = 0;
            islast = false;
        }
    }
    return data;
}

