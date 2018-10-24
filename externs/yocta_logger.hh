/**
* @fileoverview Copyright (c) 2017,2018, by Stefano Gualandi
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*/

#pragma once

#include <string>
#include <chrono>
#include <memory>

#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <fstream>
#include <sstream>

// In order to use PRId64
#include <inttypes.h>

namespace yocta {
// Why Yocto? Miss Yocta Geek is the girlfriend of mister Yocto
// (for real, please, read the following link)
// Wikipedia: https://en.wikipedia.org/wiki/Yocto-


// Support for forwarding the content of fprintf to my logger class
// From Stackoverflow at:
// https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
template<typename ... Args>
std::string fmt(const std::string& format, Args ... args) {
   size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
   std::unique_ptr<char[]> buf(new char[size]);
   snprintf(buf.get(), size, format.c_str(), args ...);
   return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}


// Read a text file and store it in a unique string
std::string readTextFile(const std::string& filename) {
   std::ifstream in(filename, std::ios::in | std::ios::binary);
   if (in) {
      std::string doc;
      in.seekg(0, std::ios::end);
      doc.resize(in.tellg());
      in.seekg(0, std::ios::beg);
      in.read(&doc[0], doc.size());
      in.close();
      return doc;
   }
   throw std::runtime_error("File not found: " + filename);
}

// Verbosity levels
enum VerbosityLevel { ERROR = 0, WARN = 1, NOTE = 2, INFO = 3, DEBUG = 4 };

class Logger {
 public:
   // Standard c'tor
   Logger() : stream(stdout), vl(VerbosityLevel::INFO) {}
   // Log to the given filename
   Logger(const std::string& filename)
      : stream(std::fopen(filename.c_str(), "w")),
        vl(VerbosityLevel::INFO)
   {}
   // Flush, close, and de c'tor
   ~Logger() {
      fflush(stream);
      if (stream != stdout)
         std::fclose(stream);
   }

   // Rule of five: Move constructor
   Logger(Logger&& o)
      : stream(o.stream), vl(o.vl)
   {}
   // Deleted constructor (this is a singleton object)
   Logger(const Logger& o) = delete;
   Logger &operator=(const Logger&) = delete;
   Logger &operator=(Logger&&) = delete;

   // Set file stream
   void setFileStream(const std::string& filename) {
      stream = std::fopen(filename.c_str(), "w");
   }
   // Set verbosity level
   void setVerbosityLevel(VerbosityLevel verbosity) {
      vl = verbosity;
   }
   // Flush out the buffer stream
   void flush() const {
      fflush(stream);
   }
   // ---- Dump message functions ---- //
   // Note this question and answer:
   // https://stackoverflow.com/questions/20588191/error-with-variadiac-template-parameter-pack-must-be-expanded

   // ERROR
   void error(const std::string& format) const {
      dump("[ERROR]", format);
      // Flush right away for errors
      fflush(stream);
   }

   template<typename ... Args>
   void error(const std::string& format, Args ... args) const {
      dump("[ERROR]", fmt(format, std::forward<Args>(args)...));
      // Flush right away for errors
      fflush(stream);
   }

   // WARNING
   void warn(const std::string& format) const {
      if (vl >= VerbosityLevel::WARN)
         dump("[WARN ]", format);
   }

   template<typename ... Args>
   void warn(const std::string& format, Args ... args) const {
      if (vl >= VerbosityLevel::WARN)
         dump("[WARN ]", fmt(format, std::forward<Args>(args)...));
   }

   // NOTE
   void note(const std::string& format) const {
      if (vl >= VerbosityLevel::NOTE)
         dump("[NOTE ]", format);
      // Flush right away for errors
      fflush(stream);
   }

   template<typename ... Args>
   void note(const std::string& format, Args ... args) const {
      if (vl >= VerbosityLevel::NOTE)
         dump("[NOTE ]", fmt(format, std::forward<Args>(args)...));
      // Flush right away for errors
      fflush(stream);
   }

   // INFO
   void info(const std::string& format) const {
      if (vl >= VerbosityLevel::INFO)
         dump("[INFO ]", format);
   }

   template<typename ... Args>
   void info(const std::string& format, Args ... args) const {
      if (vl >= VerbosityLevel::INFO)
         dump("[INFO ]", fmt(format, std::forward<Args>(args)...));
   }

   // DEBUG
   void debug(const std::string& format) const {
      if (vl >= VerbosityLevel::DEBUG)
         dump("[DEBUG]", format);
   }

   template<typename ... Args>
   void debug(const std::string& format, Args ... args) const {
      if (vl >= VerbosityLevel::DEBUG)
         dump("[DEBUG]", fmt(format, std::forward<Args>(args)...));
   }

 private:
   // Output stream
   std::FILE* stream;
   // Verbosity level
   VerbosityLevel vl;

   // Dump the message to the stream, with datetime format
   void dump(const char* msg, const std::string& message) const {
      using namespace std;
      using namespace std::chrono;
      // get current time
      auto now = system_clock::now();
      time_t now_time = system_clock::to_time_t(now);
      // get datetime
      char date[100];
      strftime(date, sizeof(date), "%Y-%m-%d %H:%M:%S", localtime(&now_time));
      // Get milliseocnds
      auto since_epoch = now.time_since_epoch();
      auto s = duration_cast<std::chrono::seconds>(since_epoch);
      since_epoch -= s;
      milliseconds milli = duration_cast<milliseconds>(since_epoch);
      // dump the string
      fprintf(stream, "%s.%.3" PRId64 " %s %s\n", date, milli.count(), msg, message.c_str());
   }
};

} // End namespace Yocto
