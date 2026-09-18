#include "app/window.hpp"
#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>
int main(int argc,char** argv) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication app(argc,argv);
    QApplication::setApplicationName("CRUCIBLE");
    QString output;
    if(app.arguments().size()==3 && app.arguments()[1]=="--smoke-test") output=app.arguments()[2];
    Window window(output);window.show();return app.exec();
}
