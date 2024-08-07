#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <string>



class UserSettings {

    public:
                                UserSettings(void);
                                UserSettings(const UserSettings& obj) = delete;
        std::string             arguments(void);
        static UserSettings&    getUserSettings(void)
                                    {
                                    static UserSettings singleInstance;
                                    return singleInstance;
                                    }
        int                     getBurnIn(void) { return burnIn; }
        int                     getChainLength(void) { check(); return chainLength; }
        int                     getMappingFrequency(void) { return mappingFrequency; }
        int                     getNumTreads(void) { check(); return numThreads; }
        std::string             getOutputFile(void) { check(); return outFile; }
        int                     getPrintFrequency(void) { check(); return printFrequency; }
        int                     getSampleFrequency(void) { check(); return sampleFrequency; }
        std::string             getTreeFile(void) { check(); return treeFile; }
        std::string             getTsvFile(void) { check(); return tsvFile; }
        void                    initializeSettings(int argc, char* argv[]);
        void                    print(void);

    private:
        bool                    check(void);
        bool                    userSettingsRead;
        std::string             treeFile;
        std::string             tsvFile;
        std::string             outFile;
        int                     burnIn;
        int                     chainLength;
        int                     printFrequency;
        int                     sampleFrequency;
        int                     mappingFrequency;
        int                     numThreads;
};

#endif
