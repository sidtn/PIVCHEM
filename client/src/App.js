import React from 'react';
import AppRouter from './AppRouter';
import 'bootstrap/dist/css/bootstrap.min.css';
import './styles/App.css';
import { BrowserRouter } from 'react-router-dom';
import AppNavbar from "./components/Navbar";


function App() {
  return (
    <div>
      <AppNavbar />
      <BrowserRouter>
        <AppRouter />
      </BrowserRouter>
    </div>
  );
}

export default App;
