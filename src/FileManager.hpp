#ifndef FileManager_hpp
#define FileManager_hpp

#include <string>
#include <vector>
#include <cstdint>

typedef std::vector<std::pair<std::string,std::string>> DirectoryList;



class FileManager {

    public:
                                    FileManager(void) = delete;
                                    FileManager(std::string dn, bool createIfMissing = false);
        
                                    // directory information
        std::string                 getDirectoryBaseName(void) const;
        std::string                 getDirectoryPath(void) const { return directoryName; }
        size_t                      numFiles(void) const { return directoryFileContents.size(); }
        
                                    // file queries
        std::vector<std::string>    filesWithExtension(const std::string& extension) const;
        bool                        containsFile(const std::string& fileName) const;
        std::string                 fullPathForFile(const std::string& fileName) const;
        
                                    // existence checks
        static bool                 directoryExists(const std::string& path);
        static bool                 fileExists(const std::string& path);
        static bool                 pathExists(const std::string& path);
        static bool                 isDirectory(const std::string& path);
        static bool                 isRegularFile(const std::string& path);
        
                                    // directory operations
        static bool                 createDirectory(const std::string& path);
        static bool                 createDirectoryRecursive(const std::string& path);
        static bool                 createDirectoryIfNeeded(const std::string& path);
        static bool                 removeDirectory(const std::string& path);
        static bool                 removeDirectoryRecursive(const std::string& path);
        
                                    // file operations
        static bool                 removeFile(const std::string& path);
        static bool                 copyFile(const std::string& source, const std::string& destination);
        static bool                 moveFile(const std::string& source, const std::string& destination);
        static bool                 renameFile(const std::string& oldPath, const std::string& newPath);
        
                                    // file information
        static std::uintmax_t       getFileSize(const std::string& path);
        static std::string          getExtension(const std::string& path);
        static std::string          getFileName(const std::string& path);
        static std::string          getParentPath(const std::string& path);
        static std::string          absolutePath(const std::string& path);
        static std::string          joinPaths(const std::string& base, const std::string& relative);
        
                                    // diagnostics
        void                        print(void) const;
        void                        refresh(void);
    
    private:
        void                        scanDirectory(void);
        std::string                 directoryName;
        DirectoryList               directoryFileContents;
};

#endif
