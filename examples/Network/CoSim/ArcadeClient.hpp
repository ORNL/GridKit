#pragma once

#include <semaphore.h>

#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <utility>

#include <sys/shm.h>

namespace GridKit
{
  namespace CoSim
  {
    /// Simulator side of the ARCADE Unix DataBroker IPC protocol.
    class ArcadeClient
    {
    public:
      using Tags = std::array<std::string, 2>;

      ArcadeClient(const Tags& inputs, const Tags& outputs, double dt, int timeout_ms)
        : stop_(open("/stop")),
          announce_(open("/msg")),
          publish_(open("/pp_sem")),
          update_(open("/up_sem")),
          timeout_(timeout_ms)
      {
        validateTags(inputs);
        validateTags(outputs);
        published_ = attach<Points>(10618);
        updated_   = attach<Points>(10619);
        metadata_  = attach<Metadata>(10620);
        initialize(*published_, outputs, 0.0);
        initialize(*updated_, inputs, arcade_initial_value);
        *metadata_ = {2, 2, dt};
      }

      ~ArcadeClient()
      {
        if (announced_)
        {
          sem_post(stop_.get());
        }
      }

      /// Announce definitions and wait for the broker.
      void start(const std::pair<double, double>& values)
      {
        setPublished(values, 0.0);
        post(announce_.get(), "Announcing ARCADE metadata");
        announced_ = true;
        waitForUpdate();
      }

      /// Publish outputs.
      void send(const std::pair<double, double>& values, double time)
      {
        setPublished(values, time);
        post(publish_.get(), "Publishing ARCADE output");
      }

      /// Wait for the broker to supply both inputs.
      std::pair<double, double> receive(double time)
      {
        waitForUpdate();
        for (const auto& point : *updated_)
        {
          if (!std::isfinite(point.value) || point.value == arcade_initial_value || point.time != time)
          {
            throw std::runtime_error(std::string("Missing, invalid, or stale ARCADE input: ") + point.name);
          }
        }
        return {(*updated_)[0].value, (*updated_)[1].value};
      }

    private:
      // IPC layout.
      struct Point
      {
        char   name[128];
        char   type[50];
        double value;
        double time;
      };

      struct Metadata
      {
        int    published;
        int    updated;
        double time_step;
      };

      /// ARCADE's unset value default.
      static constexpr double arcade_initial_value = -1.0e14;

      static void post(sem_t* semaphore, const char* operation)
      {
        if (sem_post(semaphore) == -1)
        {
          throw std::system_error(errno, std::generic_category(), operation);
        }
      }

      struct CloseSemaphore
      {
        void operator()(sem_t* handle) const
        {
          sem_close(handle);
        }
      };

      using Semaphore = std::unique_ptr<sem_t, CloseSemaphore>;

      static Semaphore open(const char* name)
      {
        auto* handle = sem_open(name, 0);
        if (handle == SEM_FAILED)
        {
          throw std::system_error(errno, std::generic_category(), std::string("Opening ARCADE semaphore ") + name);
        }
        return Semaphore(handle);
      }

      struct Detach
      {
        void operator()(const void* address) const
        {
          shmdt(address);
        }
      };

      template <typename T>
      static std::unique_ptr<T, Detach> attach(key_t key)
      {
        const int id = shmget(key, sizeof(T), IPC_CREAT | 0600);
        if (id == -1)
        {
          throw std::system_error(errno, std::generic_category(), "Opening ARCADE shared memory");
        }
        auto* address = shmat(id, nullptr, 0);
        if (address == reinterpret_cast<void*>(-1))
        {
          throw std::system_error(errno, std::generic_category(), "Attaching ARCADE shared memory");
        }
        return std::unique_ptr<T, Detach>(static_cast<T*>(address));
      }

      using Points = std::array<Point, 2>;

      static void validateTags(const Tags& tags)
      {
        for (const auto& tag : tags)
        {
          if (tag.size() >= sizeof(Point::name) || tag.find_first_of(" \t\n\r\f\v,;:") != std::string::npos)
          {
            throw std::invalid_argument("Invalid ARCADE tag: " + tag);
          }
        }
      }

      static void initialize(Points& points, const Tags& tags, double initial_value)
      {
        points = {};
        for (std::size_t i = 0; i < points.size(); ++i)
        {
          std::memcpy(points[i].name, tags[i].c_str(), tags[i].size() + 1);
          std::strcpy(points[i].type, "DOUBLE");
          points[i].value = initial_value;
        }
      }

      void setPublished(const std::pair<double, double>& values, double time)
      {
        if (!std::isfinite(values.first) || !std::isfinite(values.second))
        {
          throw std::runtime_error("Cannot send non-finite co-simulation output");
        }
        (*published_)[0].value = values.first;
        (*published_)[1].value = values.second;
        for (auto& point : *published_)
        {
          point.time = time;
        }
      }

      static bool tryWait(sem_t* semaphore, const char* operation)
      {
        if (sem_trywait(semaphore) == 0)
        {
          return true;
        }
        if (errno == EAGAIN || errno == EINTR)
        {
          return false;
        }
        throw std::system_error(errno, std::generic_category(), operation);
      }

      void waitForUpdate()
      {
        using namespace std::chrono_literals;
        const auto deadline = std::chrono::steady_clock::now() + timeout_;
        for (;;)
        {
          if (tryWait(stop_.get(), "Checking ARCADE stop signal"))
          {
            post(stop_.get(), "Restoring ARCADE stop signal");
            throw std::runtime_error("ARCADE DataBroker stopped before the exchange completed");
          }
          if (tryWait(update_.get(), "Waiting for ARCADE update"))
          {
            return;
          }
          if (std::chrono::steady_clock::now() >= deadline)
          {
            throw std::runtime_error("ARCADE update timed out");
          }
          std::this_thread::sleep_for(1ms);
        }
      }

      Semaphore                         stop_;
      Semaphore                         announce_;
      Semaphore                         publish_;
      Semaphore                         update_;
      std::chrono::milliseconds         timeout_;
      std::unique_ptr<Points, Detach>   published_;
      std::unique_ptr<Points, Detach>   updated_;
      std::unique_ptr<Metadata, Detach> metadata_;
      bool                              announced_{false};
    };
  } // namespace CoSim
} // namespace GridKit
