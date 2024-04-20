import React from 'react';
import Login from './pages/Login';
import { Navigate, Route, Routes } from "react-router-dom";
import ChemEditor from "./pages/ChemEditor";


const AppRouter = () => {
  return (
    <Routes>
      <Route path="/editor" element={<ChemEditor/>}/>
      <Route path="/login" element={<Login/>}/>
      <Route path="/*" element={<Navigate to="/login"/>}/>
    </Routes>
  );
};

export default AppRouter;