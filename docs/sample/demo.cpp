#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

/*void traverse(int arr[3]) {
    for (auto &ele : arr) // can not traverse array pointer
        ele *= 2;
}*/

void func(int i) { std::cout << "void func(int) is called. \n"; }
void func(char *i) { std::cout << "void func(char *) is called. \n"; }

class TestClass {
public:
    TestClass(const std::string &tag) : tag(tag) { 
        std::cout << "Object " << tag << " constructed. \n"; 
    }
    ~TestClass() { std::cout << "Object " << tag << " destructed. \n"; }
private:
    const std::string tag;
};

int main() {
    // 1. Type inference: auto
    std::map<int, std::string> a;
    std::map<int, std::string>::iterator it1 = a.find(0); // explicit type declaration
    auto it2 = a.find(0); // 'auto' type inference

    // 2. Range-based for-loop
    std::vector<int> v{4, 2, 5};
    for (size_t i = 0; i < 3; i++) // index-based
        v[i] *= 2;
    for (auto &ele : v) // range-based
        ele *= 2;

    int arr[3] = {3, 5, 7}; // can also traverse T[]
    for (auto &ele : arr) 
        ele -= 1;

    // 3. Smart pointer with exclusive ownership: std::unique_ptr
    {
        std::unique_ptr<TestClass> ptr(new TestClass("1"));
        ptr = std::unique_ptr<TestClass>(new TestClass("2")); // automatically delete "1"
    } // automatically delete "2" when stepping out of scope

    // 4. nullptr instead of NULL
    func(NULL);
    func(nullptr);
}