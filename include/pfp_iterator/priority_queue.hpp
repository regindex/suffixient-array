#include <iostream>
#include <queue>
#include <vector>
#include <functional>

template <typename T>
class PriorityQueue {
public:

    PriorityQueue()
        : pq(comp) {}

    void push(const T& item)
    {
        pq.push(item);
    }

    void pop()
    {
        if (!pq.empty())
            pq.pop();
    }

    T top() const {
        if (!pq.empty()) {
            return pq.top();
        }
        throw std::runtime_error("Priority queue is empty");
    }

    bool empty() const
    {
        return pq.empty();
    }

    size_t size() const
    {
        return pq.size();
    }

private:

    std::function<bool(const T&, const T&)> comp = [](const T& a, const T& b){ return *a.first > *b.first; };
    std::priority_queue<T, std::vector<T>, std::function<bool(const T&, const T&)>> pq;
};