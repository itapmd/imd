/* default return codes for buttons */
#define CANCEL_BUTTON   0
#define OK_BUTTON       1
#define ALL_BUTTON      2
#define ENTRIES(x)              (sizeof((x))/sizeof((x)[0]))

/* some useful types */
typedef struct   /* a dialogs buttons */
{
  char *Name,                  /* the widgets name */
       *Label;                 /* the buttons text */
  XtCallbackProc  Callback;    /* the buttons select-handler */
  XtPointer ClientData;        /* data passed to the handler */
  Widget W;                    /* the button widget itself */
} DialogButtonType, *DialogButtonTypePtr;

struct radio_struct {
  int  *variable;
  int  value;
};

static void CB_OK(Widget W, XtPointer ClientData, XtPointer CallData);

void StartDialog(Widget Popup);
void EndDialog(Widget Popup);
void MessageDialog(char *MessageText);
int QuestionDialog(char *MessageText);
int DialogEventLoop(int *Code);
static void CB_OK(Widget W, XtPointer ClientData, XtPointer CallData);
Widget AddButtons(Widget Parent, Widget Top,
                  DialogButtonTypePtr Buttons, size_t Count);
static void CenterDialog(Widget Popup, Boolean CenterX, Boolean CenterY);
