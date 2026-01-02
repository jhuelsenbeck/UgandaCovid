#ifndef MetaData_hpp
#define MetaData_hpp

#include <map>
#include <set>
#include <string>
#include <vector>
class Node;
class Tree;

struct CollectionDate {

    unsigned day;
    unsigned month;
    int year;
};

typedef std::map<std::string,std::vector<std::string>> name_map;
typedef std::map<std::string,int> area_map;
typedef std::map<std::string,CollectionDate> time_map;
typedef std::map<std::pair<int,int>,int> interval_map;


class MetaData {

    public:
                                    MetaData(void) = delete;
                                    MetaData(std::string fileName, std::string rd);
        void                        assignNodeTimes(Tree* t);
        void                        assignTreeTipInfo(Tree* t);
        void                        assignTimeIntervals(Tree* t, std::vector<std::string> boundaryDates);
        int                         daysFromCivil(int y, unsigned m, unsigned d);
        std::vector<std::string>    getAreas(void);
        area_map&                   getAreaMap(void) { return areas; }
        time_map&                   getCollectionDates(void) { return collectionDates; }
        int                         getIntervalId(double t);
        interval_map&               getIntervalMap(void) { return intervalInfo; }
        int                         getNumAreas(void) { return (int)areas.size(); }
        CollectionDate              getRootDate(void) { return rootDate; }
        void                        print(void);
    
    private:
        std::string                 extractAreaInfo(std::vector<std::string>& vec);
        CollectionDate              extractDateInfo(std::string str);
        void                        incrementIntervalTimes(Node* p, std::vector<double>& intervalDurations);
        double                      iterateBranchTimes(Tree* t);
        double                      iterateBranchTimesUp(Tree* t);
        double                      pickBestTime(Node* p);
        void                        removeMissingAreas(Tree* t);
        double                      sumSquares(Tree* t);
        void                        tipToRootInfo(Tree* t, std::string fn);
        void                        tokenizeString(std::string& str, std::vector<std::string>& tokens, std::string& key);
        name_map                    values;
        area_map                    areas;
        time_map                    collectionDates;
        interval_map                intervalInfo;
        CollectionDate              earliestDate;
        CollectionDate              rootDate;
};

#endif
