
// Windows Graphical User Interface
#include <tchar.h>
#using <mscorlib.dll>
#using <System.DLL>
#using <System.Windows.Forms.DLL>
#using <System.Drawing.DLL>
#pragma comment( linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"" )
using namespace System;
using namespace System::Windows::Forms;
using namespace System::Drawing;

// Class for button control
__gc class Handlers {
public:
  Form *form;
  void button_Click(System::Object* s, System::EventArgs* e)
  {
    DialogResult res = MessageBox::Show("Exit?", "My Program",
      MessageBoxButtons::OKCancel);
    if (res == DialogResult::OK) {
      form->Close();
    }
  }
};

double mainline(double);
#ifdef _UNICODE
  int wmain(void)
#else
  int main(void)
#endif

{
  // Definitions of form objects
  double result;
  Form *myform = new Form();
  Handlers *handler = new Handlers();
  Label *lab1 = new Label();
  handler->form = myform;
  result = mainline(3.0);

  // Button definitions
  Button *b1 = new Button();
  b1->Text = "Click";
  b1->Left = 50;
  b1->Top = 30;
  b1->Click += new System::EventHandler(handler, 
   Handlers::button_Click);

  // Form definitions
  myform->Width = 500;
  myform->Height = 400;
  myform->Text = "My Form";

  // Label definitions
  lab1->AutoSize = true;
  lab1->Text = result.ToString();

  // Sequence of event handling
  myform->Controls->Add(lab1);
  myform->ShowDialog();
  return 0;
}
