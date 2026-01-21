#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "FileManager.hpp"



FileManager::FileManager(std::string dn, bool createIfMissing) : directoryName(std::move(dn)) {

    if (createIfMissing && !directoryExists(directoryName))
        {
        if (!createDirectoryRecursive(directoryName))
            throw std::runtime_error("FileManager: Failed to create directory '" + directoryName + "'");
        }
    
    if (!directoryExists(directoryName))
        throw std::runtime_error("FileManager: Directory does not exist '" + directoryName + "'");
    
    scanDirectory();
}

void FileManager::scanDirectory(void) {

    directoryFileContents.clear();
    
    for (const auto& entry : std::filesystem::directory_iterator(directoryName))
        {
        std::filesystem::path fp = entry.path();
        std::string filePath = fp.string();
        std::string fileExtension = fp.extension().string();
        
        // Remove leading dot from extension
        if (!fileExtension.empty() && fileExtension[0] == '.')
            fileExtension.erase(0, 1);
        
        directoryFileContents.emplace_back(filePath, fileExtension);
        }
    std::sort(directoryFileContents.begin(), directoryFileContents.end());
}

void FileManager::refresh(void) {

    scanDirectory();
}

std::string FileManager::getDirectoryBaseName(void) const {

    std::filesystem::path p(directoryName);
    if (p.filename().empty())
        p = p.parent_path();
    return p.filename().string();
}

std::vector<std::string> FileManager::filesWithExtension(const std::string& extension) const {

    std::vector<std::string> fileList;
    for (const auto& f : directoryFileContents)
        {
        if (f.second == extension)
            fileList.push_back(f.first);
        }
    return fileList;
}

bool FileManager::containsFile(const std::string& fileName) const {

    for (const auto& f : directoryFileContents)
        {
        std::filesystem::path p(f.first);
        if (p.filename().string() == fileName)
            return true;
        }
    return false;
}

std::string FileManager::fullPathForFile(const std::string& fileName) const {

    for (const auto& f : directoryFileContents)
        {
        std::filesystem::path p(f.first);
        if (p.filename().string() == fileName)
            return f.first;
        }
    return "";
}

bool FileManager::directoryExists(const std::string& path) {

    return std::filesystem::exists(path) && std::filesystem::is_directory(path);
}

bool FileManager::fileExists(const std::string& path) {

    return std::filesystem::exists(path) && std::filesystem::is_regular_file(path);
}

bool FileManager::pathExists(const std::string& path) {

    return std::filesystem::exists(path);
}

bool FileManager::isDirectory(const std::string& path) {

    return std::filesystem::is_directory(path);
}

bool FileManager::isRegularFile(const std::string& path) {

    return std::filesystem::is_regular_file(path);
}

bool FileManager::createDirectory(const std::string& path) {

    try
        {
        return std::filesystem::create_directory(path);
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::createDirectoryRecursive(const std::string& path) {

    try
        {
        return std::filesystem::create_directories(path);
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::createDirectoryIfNeeded(const std::string& path) {

    if (directoryExists(path))
        return true;
    return createDirectoryRecursive(path);
}

bool FileManager::removeDirectory(const std::string& path) {

    try
        {
        if (!directoryExists(path))
            return false;
        return std::filesystem::remove(path);
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::removeDirectoryRecursive(const std::string& path) {

    try
        {
        if (!directoryExists(path))
            return false;
        std::filesystem::remove_all(path);
        return true;
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::removeFile(const std::string& path) {

    try
        {
        if (!fileExists(path))
            return false;
        return std::filesystem::remove(path);
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::copyFile(const std::string& source, const std::string& destination) {

    try
        {
        std::filesystem::copy_file(source, destination, 
                                   std::filesystem::copy_options::overwrite_existing);
        return true;
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::moveFile(const std::string& source, const std::string& destination) {

    try
        {
        std::filesystem::rename(source, destination);
        return true;
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return false;
        }
}

bool FileManager::renameFile(const std::string& oldPath, const std::string& newPath) {

    return moveFile(oldPath, newPath);
}

std::uintmax_t FileManager::getFileSize(const std::string& path) {

    try
        {
        return std::filesystem::file_size(path);
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return static_cast<std::uintmax_t>(-1);
        }
}

std::string FileManager::getExtension(const std::string& path) {

    std::filesystem::path p(path);
    std::string ext = p.extension().string();
    if (!ext.empty() && ext[0] == '.')
        ext.erase(0, 1);
    return ext;
}

std::string FileManager::getFileName(const std::string& path) {

    return std::filesystem::path(path).filename().string();
}

std::string FileManager::getParentPath(const std::string& path) {

    return std::filesystem::path(path).parent_path().string();
}

std::string FileManager::absolutePath(const std::string& path) {

    try
        {
        return std::filesystem::absolute(path).string();
        }
    catch (const std::filesystem::filesystem_error&)
        {
        return path;
        }
}

std::string FileManager::joinPaths(const std::string& base, const std::string& relative) {

    std::filesystem::path p(base);
    p /= relative;
    return p.string();
}

void FileManager::print(void) const {

    size_t len = 0;
    for (const auto& f : directoryFileContents)
        {
        if (f.first.length() > len)
            len = f.first.length();
        }

    int n = 0;
    for (const auto& f : directoryFileContents)
        {
        std::cout << std::setw(4) << ++n << " -- " << f.first << " ";
        for (size_t i = 0; i < len - f.first.length(); i++)
            std::cout << " ";
        std::cout << f.second << std::endl;
        }
}
