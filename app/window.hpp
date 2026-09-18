#pragma once
#include <QMainWindow>
#include <memory>
class Window : public QMainWindow {
public:
    explicit Window(QString smokeOutput = {});
    ~Window();
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
