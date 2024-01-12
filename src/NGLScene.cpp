#include <QMouseEvent>
#include <QGuiApplication>
#include <QGridLayout>
#include <QPushButton>
#include <QRect>


#include "NGLScene.h"
#include <ngl/NGLInit.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include <ngl/Util.h>
#include <iostream>

NGLScene::NGLScene()
{
    // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
    setTitle("Fluid Simulator");
    m_widget = QWidget::createWindowContainer(this);
    m_widget->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);    

}


NGLScene::~NGLScene()
{
    std::cout<<"Shutting down NGL, removing VAOs and Shaders\n";
}



void NGLScene::resizeGL(int _w , int _h)
{
    m_win.width  = static_cast<int>( _w * devicePixelRatio() );
    m_win.height = static_cast<int>( _h * devicePixelRatio() );
    m_project=ngl::perspective(45.0f,static_cast<float>(_w)/_h,0.1f,100.0f);
}


void NGLScene::initializeGL()
{
    // we must call that first before any other GL commands to load and link the
    // gl commands from the lib, if that is not done program will crash
    ngl::NGLInit::initialize();
    glClearColor(0.7f, 0.7f, 0.7f, 1.0f);			   // Grey Background
    // enable depth testing for drawing
    glEnable(GL_DEPTH_TEST);
    // enable multisampling for smoother drawing
    glEnable(GL_MULTISAMPLE);

    m_emitter=std::make_unique<Simulator>(ngl::Vec3(0.0f, 10.0f, 0.0f), 10000);

    ngl::ShaderLib::loadShader("ParticleShader","shaders/ParticleVertex.glsl","shaders/ParticleFragment.glsl");
    ngl::ShaderLib::use("ParticleShader");
    m_view=ngl::lookAt({0,0,24},{0,0,0},{0,1,0});
    ngl::ShaderLib::setUniform("MVP",m_project*m_view);
    startTimer(10);
    ngl::VAOPrimitives::createLineGrid("floor",40,40,10);

    m_previousTime = std::chrono::steady_clock::now();

    QGridLayout *m_grid = new QGridLayout(m_widget);
    m_widget->setLayout(m_grid); // Example layout

    QPushButton* m_button = new QPushButton("push", m_widget);
    m_grid->addWidget(m_button);
    // set the name
    m_button->setObjectName(QString::fromUtf8("button"));
    // set the geometry
    m_button->setGeometry(QRect(10, 80, 100, 32));
    // set the text of the button
    m_button->setText("push");
}



void NGLScene::paintGL()
{
    // clear the screen and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0,0,m_win.width,m_win.height);
    auto rotX = ngl::Mat4::rotateX(m_win.spinXFace);
    auto rotY = ngl::Mat4::rotateY(m_win.spinYFace);
    auto mouseRotation= rotX* rotY;
    mouseRotation.m_m[3][0] = m_modelPos.m_x;
    mouseRotation.m_m[3][1] = m_modelPos.m_y;
    mouseRotation.m_m[3][2] = m_modelPos.m_z;


    ngl::ShaderLib::use("ParticleShader");
    ngl::ShaderLib::setUniform("MVP",m_project*m_view*mouseRotation);
    m_emitter->draw();
    ngl::ShaderLib::use(ngl::nglColourShader);
    ngl::ShaderLib::setUniform("Colour",1.0f,0.0f,1.0f,1.0f);
    ngl::ShaderLib::setUniform("MVP",m_project*m_view*mouseRotation);
    ngl::VAOPrimitives::draw("floor");


//  ngl::VAOPrimitives::draw("teapot");

}

//----------------------------------------------------------------------------------------------------------------------

void NGLScene::keyPressEvent(QKeyEvent *_event)
{
    // this method is called every time the main window recives a key event.
    // we then switch on the key value and set the camera in the GLWindow
    switch (_event->key())
    {
        // escape key to quite
        case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
        case Qt::Key_Space :
            m_win.spinXFace=0;
            m_win.spinYFace=0;
            m_modelPos.set(ngl::Vec3::zero());

            break;
        default : break;
    }
    // finally update the GLWindow and re-draw

    update();
}

void NGLScene::timerEvent(QTimerEvent *_event)
{
    auto now = std::chrono::steady_clock::now();
    auto delta = std::chrono::duration<float,std::chrono::seconds::period>(now-m_previousTime).count();
    m_previousTime=now;
    //std::cout<<delta.count()<<"\n";
    m_emitter->update(delta);
    update();
}