#ifndef MetaData_hpp
#define MetaData_hpp

#include <map>
#include <set>
#include <string>
#include <vector>
class Tree;

typedef std::map<std::string,std::vector<std::string>> name_map;
typedef std::map<std::string,int> area_map;



class MetaData {

    public:
                                    MetaData(void) = delete;
                                    MetaData(std::string fileName);
        void                        assignTreeTipAreas(Tree* t);
        std::vector<std::string>    getAreas(void);
        int                         getNumAreas(void) { return (int)areas.size(); }
        void                        print(void);
    
    private:
        std::string                 extractAreaInfo(std::vector<std::string>& vec);
        void                        tokenizeString(std::string& str, std::vector<std::string>& tokens, std::string& key);
        name_map                    values;
        area_map                    areas;
};

#endif
