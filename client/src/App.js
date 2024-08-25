import "ketcher-react/dist/index.css";
import "./App.css";
import 'react-toastify/dist/ReactToastify.css';

import React from "react";
import AppNavbar from "./components/AppNavbar";
import { BrowserRouter } from "react-router-dom";
import { AppRouter } from "./components/AppRouter";
import { ToastContainer } from "react-toastify";

function App() {
  return (
    <div className="app">
      <BrowserRouter>
        <AppNavbar />
        <AppRouter />
      </BrowserRouter>
      <ToastContainer />
    </div>
  );
}

export default App;
