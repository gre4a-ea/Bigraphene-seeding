#include <iostream>
#include <vector>
#include <string>


int main()
{
    std::vector <int> vect;
    std::string string;

    getline(std::cin, string);
    string += ' ';

    while (string.length() > 0){
        std::string tmp = string.substr(0, string.find(' '));
        string.erase(0, string.find(' ') + 1);
        int tmp_int = stoi(tmp);
        vect.push_back(tmp_int);
      }

    int n = vect.size();
    int i = 0;
    while (i < n){
        int j = i;
        while (j > 0 && vect[j - 1] > vect[j])
        {
            std::swap(vect[j], vect[j - 1]);
            j -= 1;
        }

        i += 1;
    }

    std::cout << "\nОтсортированный массив:\n";

    for (int i = 0; i < n; i++)
        std::cout << vect[i] << ' ';

    std::cout << '\n';

    return 0;
}
