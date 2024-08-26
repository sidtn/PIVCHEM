import { toast } from "react-toastify";


export default function notify(message, action='error') {
  
  let style = {
    position: "bottom-right",
    autoClose: 3000,
    hideProgressBar: false,
    closeOnClick: true,
    pauseOnHover: true,
    draggable: true,
    progress: undefined,
    theme: "light",
  }

  if (action === 'error') {
    toast.error(message, style)
  } else {
    toast.success(message, style)
  }

}
