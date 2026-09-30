#include "app/window.hpp"
#include "core/session.hpp"
#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QFile>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QTimer>
#include <QVBoxLayout>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <cmath>
#include <algorithm>

struct Window::Impl {
    Window* window;
    crucible::Definition definition;
    std::unique_ptr<crucible::Session> session;
    QVTKOpenGLNativeWidget* view;
    QLabel *status,*metrics;
    QDoubleSpinBox* pressure;
    QComboBox* field;
    vtkNew<vtkGenericOpenGLRenderWindow> render;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkPolyData> mesh;
    vtkNew<vtkDoubleArray> values;
    vtkNew<vtkPolyDataMapper> mapper;
    vtkNew<vtkLookupTable> colors;
    vtkNew<vtkActor> actor;
    vtkNew<vtkScalarBarActor> legend;
    std::uint64_t generation=0,control=0,requestedControl=0;
    QString smoke;
    QElapsedTimer elapsed;
    int phase=0;
    explicit Impl(Window* w,QString output):window(w),smoke(output) {
        auto root=new QWidget;auto layout=new QVBoxLayout(root);
        auto title=new QLabel("CRUCIBLE  /  Gas-flow laboratory");
        title->setStyleSheet("font-size: 23px; font-weight: 600; padding: 8px;");layout->addWidget(title);
        auto note=new QLabel("Axisymmetric nozzle · ideal gas · inviscid slip walls · prepared flowing initial state\nGas-dynamics prototype: no combustion, turbulence, or real-engine validation yet.");
        note->setStyleSheet("padding: 8px;");layout->addWidget(note);
        auto controls=new QHBoxLayout;
        auto run=new QPushButton("Run");auto pause=new QPushButton("Pause");auto reset=new QPushButton("Restart");
        run->setObjectName("run");pause->setObjectName("pause");
        controls->addWidget(run);controls->addWidget(pause);controls->addWidget(reset);
        controls->addSpacing(20);controls->addWidget(new QLabel("Inlet reservoir pressure"));
        pressure=new QDoubleSpinBox;pressure->setRange(240,360);pressure->setValue(300);pressure->setSuffix(" kPa");pressure->setSingleStep(10);
        controls->addWidget(pressure);auto apply=new QPushButton("Apply live");apply->setObjectName("applyPressure");controls->addWidget(apply);controls->addStretch();
        field=new QComboBox;field->addItems({"Mach number","Pressure (kPa)","Temperature (K)","Density (kg/m³)"});controls->addWidget(field);
        auto fit=new QPushButton("Fit view");controls->addWidget(fit);
        QObject::connect(fit,&QPushButton::clicked,w,[this]{renderer->ResetCamera();renderer->GetActiveCamera()->SetParallelScale(.15);render->Render();});
        layout->addLayout(controls);
        view=new QVTKOpenGLNativeWidget;view->setMinimumHeight(350);view->setRenderWindow(render);render->AddRenderer(renderer);
        renderer->SetBackground(.065,.09,.12);
        vtkNew<vtkInteractorStyleImage> navigation;
        render->GetInteractor()->SetInteractorStyle(navigation);
        layout->addWidget(view,1);
        auto caption=new QLabel("Axial distance 0–0.60 m · radius ±0.035 m · mirrored axisymmetric section · true geometric proportions");layout->addWidget(caption);
        metrics=new QLabel;metrics->setStyleSheet("font-size: 15px; padding: 14px; background: #eef2f5; color: #172636;");layout->addWidget(metrics);
        status=new QLabel;layout->addWidget(status);auto save=new QPushButton("Save current field as CSV…");layout->addWidget(save,0,Qt::AlignRight);
        window->setCentralWidget(root);window->resize(1140,700);window->setWindowTitle("CRUCIBLE — nozzle experiment");
        colors->SetNumberOfTableValues(256);colors->Build();
        const double anchors[5][3]={{.267,.005,.329},{.230,.322,.546},{.128,.567,.551},{.369,.789,.383},{.993,.906,.144}};
        for(int i=0;i<256;++i){double x=i*4.0/255;int j=std::min(3,static_cast<int>(x));double f=x-j;
            colors->SetTableValue(i,anchors[j][0]*(1-f)+anchors[j+1][0]*f,anchors[j][1]*(1-f)+anchors[j+1][1]*f,anchors[j][2]*(1-f)+anchors[j+1][2]*f); }
        mapper->SetInputData(mesh);mapper->SetLookupTable(colors);mapper->SetScalarModeToUseCellData();actor->SetMapper(mapper);renderer->AddActor(actor);
        legend->SetLookupTable(colors);legend->SetTitle("Mach");legend->SetNumberOfLabels(5);legend->SetWidth(.09);legend->SetHeight(.65);legend->SetPosition(.89,.18);renderer->AddActor(legend);
        QObject::connect(run,&QPushButton::clicked,w,[this]{session->run();});
        QObject::connect(pause,&QPushButton::clicked,w,[this]{session->pause();});
        QObject::connect(reset,&QPushButton::clicked,w,[this]{restart();});
        QObject::connect(apply,&QPushButton::clicked,w,[this]{requestedControl=session->setTotalPressure(pressure->value()*1000);});
        QObject::connect(field,&QComboBox::currentIndexChanged,w,[this]{generation=0;tick();});
        QObject::connect(save,&QPushButton::clicked,w,[this]{
            auto snapshot=session->latest();if(!snapshot)return;
            auto path=QFileDialog::getSaveFileName(window,"Save field",{},"CSV (*.csv)");if(path.isEmpty())return;
            QFile file(path);if(!file.open(QIODevice::WriteOnly)){status->setText("Could not save field.");return;}
            QByteArray data="time_s,z_m,r_m,density_kg_m3,axial_velocity_m_s,radial_velocity_m_s,pressure_Pa\n";
            crucible::Mesh geometry(snapshot->definition);
            for(std::size_t k=0;k<snapshot->cells.size();++k){auto q=snapshot->cells[k];auto c=geometry.cells[k];
                data+=QString("%1,%2,%3,%4,%5,%6,%7\n").arg(snapshot->measurements.time,0,'g',16).arg(c.z,0,'g',16).arg(c.r,0,'g',16).arg(q.rho,0,'g',16).arg(q.uz,0,'g',16).arg(q.ur,0,'g',16).arg(q.p,0,'g',16).toUtf8();}
            if(file.write(data)!=data.size())status->setText("Field save failed.");
        });
        restart();auto timer=new QTimer(w);QObject::connect(timer,&QTimer::timeout,w,[this]{tick();});timer->start(50);elapsed.start();
    }
    void restart(){session=std::make_unique<crucible::Session>(definition);generation=0;requestedControl=0;pressure->setValue(300);}
    void geometry(const crucible::FieldSnapshot& s){
        vtkNew<vtkPoints> points;vtkNew<vtkCellArray> faces;
        const auto& d=s.definition;
        for(int sign:{-1,1})for(int i=0;i<d.nz;++i)for(int j=0;j<d.nr;++j){
            vtkIdType ids[4];int k=0;
            for(auto corner:{std::pair{i,j},std::pair{i+1,j},std::pair{i+1,j+1},std::pair{i,j+1}})
                ids[k++]=points->InsertNextPoint(d.length*corner.first/d.nz,sign*s.radius[corner.first]*corner.second/d.nr,0);
            faces->InsertNextCell(4,ids);
        }
        mesh->SetPoints(points);mesh->SetPolys(faces);values->SetNumberOfTuples(2*s.cells.size());mesh->GetCellData()->SetScalars(values);
        renderer->GetActiveCamera()->ParallelProjectionOn();renderer->ResetCamera();renderer->GetActiveCamera()->SetParallelScale(.15);
    }
    void tick(){
        if(!smoke.isEmpty() && elapsed.elapsed()>20000){QApplication::exit(2);return;}
        auto s=session->latest();if(!s){
            if(session->status()==crucible::RunState::Failed)status->setText(QString::fromStdString(session->error()));
            return;
        }
        if(!mesh->GetNumberOfPoints())geometry(*s);
        auto m=s->measurements;
        if(s->generation!=generation){generation=s->generation;double lo=1e300,hi=-1e300;
            for(std::size_t k=0;k<s->cells.size();++k){auto q=s->cells[k];double v=0;
                switch(field->currentIndex()){case 0:v=std::hypot(q.uz,q.ur)/std::sqrt(definition.gas.gamma*q.p/q.rho);break;case 1:v=q.p/1000;break;case 2:v=q.p/(q.rho*definition.gas.specificR);break;default:v=q.rho;}
                values->SetValue(k,v);values->SetValue(k+s->cells.size(),v);lo=std::min(lo,v);hi=std::max(hi,v);}
            if(field->currentIndex()==0){lo=0;hi=3;}else if(hi-lo<1e-10)hi=lo+1;
            mapper->SetScalarRange(lo,hi);colors->SetRange(lo,hi);legend->SetTitle(field->currentText().toUtf8().constData());values->Modified();render->Render();
            metrics->setText(QString("Physical time  %1 ms     |     Inlet  %2 kg/s     Outlet  %3 kg/s     |     Exit Mach  %4\n"
                "Device thrust  %5 N  =  supply plane %6 N  +  walls %7 N  −  ambient %8 N     (exit-plane estimate %9 N)\n"
                "Mass / energy / axial-momentum balance error  %10 / %11 / %12")
                .arg(m.time*1000,0,'f',3).arg(m.inletMassFlow,0,'f',4).arg(m.outletMassFlow,0,'f',4).arg(m.exitMach,0,'f',3)
                .arg(m.deviceThrust,0,'f',2).arg(m.inletMomentumFlux,0,'f',2).arg(m.wallAxialForce,0,'f',2).arg(m.ambientAxialForce,0,'f',2)
                .arg(m.exitPlaneThrust,0,'f',2).arg(m.massBalanceError,0,'e',1).arg(m.energyBalanceError,0,'e',1).arg(m.momentumBalanceError,0,'e',1));
        }
        if(m.steps==0) metrics->setText("Prepared initial state · press Run to measure boundary flow and device thrust.");
        QString state="Initializing";switch(session->status()){case crucible::RunState::Paused:state="Paused";break;case crucible::RunState::Running:state="Running";break;case crucible::RunState::PauseRequested:state="Pausing";break;case crucible::RunState::Failed:state="Failed: "+QString::fromStdString(session->error());break;default:break;}
        status->setText(QString("%1 · %2 accepted steps · applied inlet %3 kPa · control %4").arg(state).arg(m.steps).arg(s->appliedTotalPressure/1000).arg(s->appliedControlSequence));
        if(requestedControl>s->appliedControlSequence)
            status->setText(status->text()+QString(" · pending control %1").arg(requestedControl));
        if(!smoke.isEmpty()){
            if(session->status()==crucible::RunState::Failed || elapsed.elapsed()>20000){QApplication::exit(2);return;}
            if(phase==0){window->findChild<QPushButton*>("run")->click();phase=1;}
            else if(phase==1 && m.steps>20){pressure->setValue(330);window->findChild<QPushButton*>("applyPressure")->click();control=requestedControl;phase=2;}
            else if(phase==2 && s->appliedControlSequence==control && m.time>.002){window->findChild<QPushButton*>("pause")->click();phase=3;}
            else if(phase==3 && session->status()==crucible::RunState::Paused){phase=4;QTimer::singleShot(200,window,[this]{QApplication::exit(window->grab().save(smoke)?0:3);});}
        }
    }
};
Window::Window(QString output):impl_(std::make_unique<Impl>(this,output)){}
Window::~Window()=default;
