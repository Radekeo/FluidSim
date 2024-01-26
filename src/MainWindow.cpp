#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QSlider>

MainWindow::MainWindow(QWidget *parent) :QMainWindow(parent), m_ui(new Ui::MainWindow)
{
    m_ui->setupUi(this);
    m_gl=new  NGLScene(this);
    m_ui->s_mainWindowGridLayout->addWidget(m_gl,0,0,2,1);

    // Connect the slider value changed signal to a slot
    connect(m_ui->vSlider, &QSlider::valueChanged, this, &MainWindow::onViscositySliderChanged);
}

void MainWindow::onViscositySliderChanged(int value) {
    float viscosity = value; // Map slider value to desired range
    m_gl->setViscosity(viscosity); // Call a method in your NGLScene to update viscosity
}
MainWindow::~MainWindow()
{
    delete m_ui;
}