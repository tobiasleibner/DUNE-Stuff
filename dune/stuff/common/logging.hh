/**
   *  \file logging.hh
   *  \brief  logging
   **/
#ifndef LOGGING_HH_INCLUDED
#define LOGGING_HH_INCLUDED

#include <fstream>
#include <ostream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <map>
#include <assert.h>
#include "misc.hh"
#include "filesystem.hh"
#include "logstreams.hh"

namespace Dune {
namespace Stuff {
namespace Common {
class Logging;
Logging& Logger();

/** \brief handles all logging
  **/
class Logging
{
protected:
  Logging()
    : matlabLogStreamPtr(0) {
    streamIDs_.push_back(LOG_ERR);
    streamIDs_.push_back(LOG_DEBUG);
    streamIDs_.push_back(LOG_INFO);
  }

  ~Logging() {
    IdVecCIter it = streamIDs_.end();

    for ( ; it != streamIDs_.begin(); --it)
    {
      delete streammap_[*it];
      streammap_[*it] = 0;
    }

    if ( (logflags_ & LOG_FILE) != 0 )
    {
      logfile_ << '\n' << Dune::Stuff::Common::String::fromTime() << ": LOG END" << std::endl;
      logfileWoTime_ << std::endl;
      logfile_.close();
      logfileWoTime_.close();
    }

    // delete the MatlabLogStream
    matlabLogStreamPtr->flush();
    matlabLogFile_.close();
    delete matlabLogStreamPtr;
    matlabLogStreamPtr = 0;
  }

public:
  /** \brief setup loglevel, logfilename
     *  \param logflags any OR'd combination of flags
     *  \param logfile filename for log, can contain paths, but creation will fail if dir is non-existant
     **/
  void create( unsigned int logflags = LOG_FILE | LOG_CONSOLE | LOG_ERR,
               std::string logfile = "dune_stuff_log",
               std::string datadir = "data",
               std::string logdir = std::string("") ) {
    logflags_ = logflags;

    // if we get a logdir from parameters append path seperator, otherwise leave empty
    // enables us to use logdir unconditionally further down
    if ( !datadir.empty() )
      datadir += "/";
    if ( !logdir.empty() )
      logdir += "/";
    logdir = datadir + logdir;

    filename_ = logdir + logfile + "_time.log";
    Stuff::Common::Filesystem::testCreateDirectory(logdir);     // could assert this if i figure out why errno is !=
                                                                // EEXIST
    filenameWoTime_ = logdir + logfile + ".log";
    if ( (logflags_ & LOG_FILE) != 0 )
    {
      logfile_.open( filename_.c_str() );
      assert( logfile_.is_open() );
      logfileWoTime_.open( filenameWoTime_.c_str() );
      assert( logfileWoTime_.is_open() );
    }
    IdVecCIter it = streamIDs_.begin();
    for ( ; it != streamIDs_.end(); ++it)
    {
      flagmap_[*it] = logflags;
      streammap_[*it] = new LogStream(*it, flagmap_[*it], logfile_, logfileWoTime_);
    }
    // create the MatlabLogStream
    std::string matlabLogFileName = logdir + logfile + "_matlab.m";
    Stuff::Common::Filesystem::testCreateDirectory(matlabLogFileName);           // could assert this if i figure out
                                                                                 // why errno is != EEXIST
    matlabLogFile_.open( matlabLogFileName.c_str() );
    assert( matlabLogFile_.is_open() );
    matlabLogStreamPtr = new MatlabLogStream(LOG_FILE, logflags_, matlabLogFile_);
  } // Create

  // ! \attention This will probably not do wht we want it to!
  void setPrefix(std::string prefix) {
    // / begin dtor
    IdVecCIter it = streamIDs_.end();

    for ( ; it != streamIDs_.begin(); --it)
    {
      delete streammap_[*it];
      streammap_[*it] = 0;
    }

    if ( (logflags_ & LOG_FILE) != 0 )
    {
      logfile_ << '\n' << Dune::Stuff::Common::String::fromTime() << ": LOG END" << std::endl;
      logfileWoTime_ << std::endl;
      logfile_.close();
      logfileWoTime_.close();
    }

    // delete the MatlabLogStream
    matlabLogStreamPtr->flush();
    matlabLogFile_.close();
    delete matlabLogStreamPtr;
    matlabLogStreamPtr = 0;
    // / end dtor

    create(logflags_, prefix);
  } // SetPrefix

  void setStreamFlags(LogFlags stream, int flags) {
    assert( stream & (LOG_ERR | LOG_INFO | LOG_DEBUG) );
    // this might result in logging to diff targtes, so we flush the current targets
    flagmap_[stream] = flags;
  }

  int getStreamFlags(LogFlags stream) {
    assert( flagmap_.find(stream) != flagmap_.end() );
    return flagmap_[stream];
  }

  template< typename Pointer, class Class >
  void log(void (Class::* pf)(std::ostream &) const, Class& c, LogFlags stream) {
    assert( flagmap_.find(stream) != flagmap_.end() );
    if ( (flagmap_[stream] & LOG_CONSOLE) != 0 )
      (c.*pf)(std::cout);
    if ( (flagmap_[stream] & LOG_FILE) != 0 )
      (c.*pf)(logfile_);
  } // Log

  template< class Class, typename Pointer >
  void log(Pointer pf, Class& c, LogFlags stream) {
    assert( stream & (LOG_ERR | LOG_INFO | LOG_DEBUG) );
    if ( (flagmap_[stream] & LOG_CONSOLE) != 0 )
      (c.*pf)(std::cout);
    if ( (flagmap_[stream] & LOG_FILE) != 0 )
    {
      (c.*pf)(logfile_);
      (c.*pf)(logfileWoTime_);
    }
  } // Log

  /** \}
     */

  /** \name Log funcs for basic types/classes
     * \{
     */
  template< class Class >
  void log(Class c, LogFlags stream) {
    assert( stream & (LOG_ERR | LOG_INFO | LOG_DEBUG) );
    if ( (flagmap_[stream] & LOG_CONSOLE) != 0 )
      std::cout << c;
    if ( (flagmap_[stream] & LOG_FILE) != 0 )
    {
      logfile_ << c;
      logfileWoTime_ << c;
    }
  } // Log

  /** \}
     */

  LogStream& getStream(int stream) {
    assert(streammap_[(LogFlags) stream]);
    return *streammap_[(LogFlags) stream];
  }

  LogStream& err() { return getStream(LOG_ERR); }
  LogStream& info() { return getStream(LOG_INFO); }
  LogStream& dbg() { return getStream(LOG_DEBUG); }
  MatlabLogStream& matlab() { return *matlabLogStreamPtr; }

  void flush() {
    for (StreamMap::iterator it = streammap_.begin();
         it != streammap_.end();
         ++it)
    {
      it->second->flush();
    }
  } // Flush

  int addStream(int flags) {
// assert( streamIDs_.find( streamID ) == streamIDs_.end() );
    static int streamID_int = 16;

    streamID_int <<= 2;
    LogFlags streamID = (LogFlags) streamID_int;
    streamIDs_.push_back(streamID);
    flagmap_[streamID] = flags | streamID;
    streammap_[streamID] = new LogStream(streamID, flagmap_[streamID], logfile_, logfileWoTime_);
    return streamID_int;
  } // AddStream

  void resume(LogStream::PriorityType prio = LogStream::default_suspend_priority) {
    for (StreamMap::iterator it = streammap_.begin();
         it != streammap_.end();
         ++it)
    {
      it->second->resume(prio);
    }
  } // Resume

  void suspend(LogStream::PriorityType prio = LogStream::default_suspend_priority) {
    for (StreamMap::iterator it = streammap_.begin();
         it != streammap_.end();
         ++it)
    {
      it->second->suspend(prio);
    }
  } // Suspend

  struct SuspendLocal
  {
    LogStream::PriorityType prio_;
    SuspendLocal(LogStream::PriorityType prio = LogStream::default_suspend_priority)
      : prio_(prio) {
      Logger().suspend(prio_);
    }

    ~SuspendLocal() {
      Logger().resume(prio_);
    }
  };

  struct ResumeLocal
  {
    LogStream::PriorityType prio_;
    ResumeLocal(LogStream::PriorityType prio = LogStream::default_suspend_priority)
      : prio_(prio) {
      Logger().resume(prio_);
    }

    ~ResumeLocal() {
      Logger().suspend(prio_);
    }
  };

private:
  std::string filename_;
  std::string filenameWoTime_;
  std::ofstream logfile_;
  std::ofstream logfileWoTime_;
  std::ofstream matlabLogFile_;
  typedef std::map< LogFlags, int > FlagMap;
  FlagMap flagmap_;
  typedef std::map< LogFlags, LogStream* > StreamMap;
  StreamMap streammap_;
  typedef std::vector< LogFlags >                 IdVec;
  typedef std::vector< LogFlags >::const_iterator IdVecCIter;
  IdVec streamIDs_;
  int logflags_;
  MatlabLogStream* matlabLogStreamPtr;

  friend Logging& Logger();
  // satisfy stricter warnings wrt copying
  Logging(const Logging&);
  Logging& operator=(const Logging&);
};

// !global Logging instance
Logging& Logger() {
  static Logging log;

  return log;
}
} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // ifndef LOGGING_HH_INCLUDED

/** Copyright (c) 2012, Rene Milk
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above copyright notice, this
   *    list of conditions and the following disclaimer.
   * 2. Redistributions in binary form must reproduce the above copyright notice,
   *    this list of conditions and the following disclaimer in the documentation
   *    and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   * The views and conclusions contained in the software and documentation are those
   * of the authors and should not be interpreted as representing official policies,
   * either expressed or implied, of the FreeBSD Project.
   **/
