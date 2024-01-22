#ifndef NGLSCENE_H_
#define NGLSCENE_H_

#include <ngl/Vec3.h>
#include "WindowParams.h"
#include <QEvent>
#include <QResizeEvent>
#include <QOpenGLWidget>
#include <memory>

#include "Simulator.h"
#include <ngl/Mat4.h>


class NGLScene : public QOpenGLWidget
{
Q_OBJECT        // must include this if you use Qt signals/slots
public :
    /// @brief Constructor for GLWindow
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Constructor for GLWindow
    /// @param [in] _parent the parent window to create the GL context in
    //----------------------------------------------------------------------------------------------------------------------
    NGLScene(QWidget *_parent );

    /// @brief dtor
    ~NGLScene() override;
public slots :

private :

protected:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the windows params such as mouse and rotations etc
    //----------------------------------------------------------------------------------------------------------------------
    WinParams m_win;
    /// @brief  The following methods must be implemented in the sub class
    /// this is called when the window is created
    void initializeGL() override;

    /// @brief this is called whenever the window is re-sized
    /// @param[in] _w the width of the resized window
    /// @param[in] _h the height of the resized window
    void resizeGL(int _w , int _h) override;
    /// @brief this is the main gl drawing routine which is called whenever the window needs to
    // be re-drawn
    void paintGL() override;

    /// @brief our model position
    ngl::Vec3 m_modelPos;
    std::unique_ptr<Simulator> m_emitter;
    ngl::Mat4 m_view;
    ngl::Mat4 m_project;
    std::chrono::steady_clock::time_point m_previousTime;

private :
    /// @brief Qt Event called when a key is pressed
    /// @param [in] _event the Qt event to query for size etc
    //----------------------------------------------------------------------------------------------------------------------
    void keyPressEvent(QKeyEvent *_event) override;
    /// @brief this method is called every time a mouse is moved
    /// @param _event the Qt Event structure

    void mouseMoveEvent (QMouseEvent * _event   ) override;
    /// @brief this method is called everytime the mouse button is pressed
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure

    void mousePressEvent ( QMouseEvent *_event  ) override;

    /// @brief this method is called everytime the mouse button is released
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    void mouseReleaseEvent (QMouseEvent *_event ) override;
    void wheelEvent( QWheelEvent* _event ) override;
    void timerEvent(QTimerEvent *_event) override;



};

#endif