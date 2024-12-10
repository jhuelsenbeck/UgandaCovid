#ifndef McmcInfo_hpp
#define McmcInfo_hpp

#include <map>
#include <string>

struct AcceptInfo {

    int     numTries;
    int     numAccepts;
};

class McmcInfo {

    public:
        void                                accept(std::string& moveStr);
        void                                print(void);
        void                                reject(std::string& moveStr);
    
    private:
        std::map<std::string,AcceptInfo>    info;
};

#endif
