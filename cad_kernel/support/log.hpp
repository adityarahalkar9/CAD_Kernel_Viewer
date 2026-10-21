// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : support/Log.hpp
//  Purpose : Minimal leveled logger (header-only, no dependencies)
//
//  * Threshold filter (default Info), swappable sink
//  * Message formatting is skipped entirely when below threshold
//  * Thread-safe: a sink is never invoked concurrently
//  * The viewer installs a custom sink later (route into its debug console);
//    the kernel itself never touches windowing or GPU APIs
// ============================================================================
#pragma once

#include <cstdio>
#include <cstdint>
#include <format>
#include <functional>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>

namespace cad {

    enum class LogLevel : std::uint8_t {
        Trace = 0,
        Debug = 1,
        Info = 2,
        Warn = 3,
        Error = 4,
        Off = 5,   // threshold only: silences everything
    };

    struct LogRecord {
        LogLevel    level{LogLevel::Info};
        std::string message{};
    };

    namespace detail {

        inline constexpr std::string_view toStringView(LogLevel level) noexcept {
            switch(level) {
                case LogLevel::Trace: return "TRACE";
                case LogLevel::Debug: return "DEBUG";
                case LogLevel::Info:  return "INFO ";
                case LogLevel::Warn:  return "WARN ";
                case LogLevel::Error: return "ERROR";
                default:              return "OFF  ";
            }
        }

        // Default sink: "[LEVEL] message" to stderr.
        inline void defaultLogSink(const LogRecord& record) {
            const std::string_view tag = toStringView(record.level);
            std::fprintf(stderr, "[%.*s] %s\n",
                         static_cast<int>(tag.size()), tag.data(),
                         record.message.c_str());
        }

        // Singletons via inline-function statics: exactly one instance across every
        // translation unit that includes this header (guaranteed for inline fns).
        inline LogLevel& logThreshold() {
            static LogLevel level = LogLevel::Info;
            return level;
        }

        inline std::function<void(const LogRecord&)>& logSink() {
            static std::function<void(const LogRecord&)> sink = &defaultLogSink;
            return sink;
        }

        inline std::mutex& logMutex() {
            static std::mutex mutex;
            return mutex;
        }

    } // namespace detail

    inline void setLogLevel(LogLevel level) noexcept {
        detail::logThreshold() = level;
    }

    inline LogLevel logLevel() noexcept {
        return detail::logThreshold();
    }

    // Installs a custom sink (route into a GUI console, a file, ...).
    // Passing nullptr restores the default stderr sink. Sinks must not throw.
    inline void setLogSink(std::function<void(const LogRecord&)> sink) {
        const std::scoped_lock lock(detail::logMutex());
        detail::logSink() = sink ? std::move(sink)
            : &detail::defaultLogSink;
    }

    // Core entry point. The threshold check happens BEFORE any formatting work.
    template <typename... Args>
    inline void logMessage(LogLevel level, std::format_string<Args...> fmt, Args&&... args) {
        if(level == LogLevel::Off) {
            return;
        }
        if(static_cast<std::uint8_t>(level) <
           static_cast<std::uint8_t>(detail::logThreshold())) {
            return;
        }

        LogRecord record;
        record.level = level;
        record.message = std::format(fmt, std::forward<Args>(args)...);

        const std::scoped_lock lock(detail::logMutex());
        if(detail::logSink()) {
            detail::logSink()(record);
        }
    }

    template <typename... Args>
    inline void logTrace(std::format_string<Args...> fmt, Args&&... args) {
        logMessage(LogLevel::Trace, fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void logDebug(std::format_string<Args...> fmt, Args&&... args) {
        logMessage(LogLevel::Debug, fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void logInfo(std::format_string<Args...> fmt, Args&&... args) {
        logMessage(LogLevel::Info, fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void logWarn(std::format_string<Args...> fmt, Args&&... args) {
        logMessage(LogLevel::Warn, fmt, std::forward<Args>(args)...);
    }

    template <typename... Args>
    inline void logError(std::format_string<Args...> fmt, Args&&... args) {
        logMessage(LogLevel::Error, fmt, std::forward<Args>(args)...);
    }

} // namespace cad